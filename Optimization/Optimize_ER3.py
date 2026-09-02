import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
import random
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# =====================================================================
# 1. Top-Level Worker Functions 
# =====================================================================
def evaluate_partition(partition_dict, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_elc):
    metrics = []
    total_cnr_sum = 0
    
    # 5 Abstract Groups for Methodological Parity
    for group_id in range(5):
        group_elc_ids = [e_id for e_id, g in partition_dict.items() if g == group_id]
        group_elc_area = overlap_matrix[group_elc_ids].sum().sum() if group_elc_ids else 0
        outside_total_area = total_valid_pixels - group_elc_area
        
        cluster_benefits = []
        
        # Inner Loop: Optimize Clusters (Case-Normalized Rate)
        if group_elc_ids:
            for cluster_id in overlap_matrix.index:
                overlap_area = overlap_matrix.loc[cluster_id, group_elc_ids].sum()
                outside_spill = cluster_areas.get(cluster_id, 0) - overlap_area
                
                tpr_gain = overlap_area / group_elc_area if group_elc_area > 0 else 0
                tnr_loss = outside_spill / outside_total_area if outside_total_area > 0 else 0
                
                benefit = tpr_gain - tnr_loss
                cluster_benefits.append((cluster_id, benefit))
                    
        # Sort all clusters by highest case-normalized benefit
        cluster_benefits.sort(key=lambda x: x[1], reverse=True)
        
        # ENFORCE CONSTRAINT: Min 2 / Max 12 clusters per group
        qualifying_clusters = [c[0] for c in cluster_benefits if c[1] > 0]
        
        if len(qualifying_clusters) < 2:
            assigned_clusters = [c[0] for c in cluster_benefits[:2]]
        else:
            assigned_clusters = qualifying_clusters[:12]
        
        cluster_area = cluster_areas[cluster_areas.index.isin(assigned_clusters)].sum() if assigned_clusters else 0
        true_positives = overlap_matrix.loc[overlap_matrix.index.isin(assigned_clusters), group_elc_ids].sum().sum() if assigned_clusters else 0
        
        false_positives = cluster_area - true_positives
        false_negatives = group_elc_area - true_positives
        true_negatives = outside_total_area - false_positives
        
        tpr = true_positives / group_elc_area if group_elc_area > 0 else 0
        tnr = true_negatives / outside_total_area if outside_total_area > 0 else 0
        
        cnr = (tpr + tnr) / 2
        total_cnr_sum += cnr
        accuracy = (true_positives + true_negatives) / total_valid_pixels if total_valid_pixels > 0 else 0
        
        temp_orders = sorted([orig_to_temp_map[c] for c in assigned_clusters])
        original_elc_strings = sorted([int_to_elc[e_id] for e_id in group_elc_ids])
        
        metrics.append({
            'Abstract Group': f"Group {group_id + 1}",
            'Optimized Clusters': str(temp_orders),
            'Optimized ER3s': str(original_elc_strings),
            'Prediction Accuracy': round(accuracy, 4),
            'Case-Normalized Rate': round(cnr, 4)
        })
        
    mean_cnr = total_cnr_sum / 5
    return mean_cnr, metrics

def run_hill_climb(args):
    seed, iterations, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, elc_list, int_to_elc, elc_areas, min_area = args
    random.seed(seed)

    # FORCE 100% COVERAGE, MIN 2 ER3s, AND MINIMUM 10% AREA PER GROUP
    valid_start = False
    while not valid_start:
        current_partition = {}
        group_counts = {g: 0 for g in range(5)}
        group_areas = {g: 0 for g in range(5)}
        available_elcs = list(elc_list)
        random.shuffle(available_elcs)
        
        # Step 1: Guarantee every group gets exactly 2 ER3s to start
        for g in range(5):
            for _ in range(2):
                if available_elcs:
                    elc = available_elcs.pop()
                    current_partition[elc] = g
                    group_counts[g] += 1
                    group_areas[g] += elc_areas[elc]

        # Step 2: Distribute the remainder, enforcing the Max 12 limit
        for elc in available_elcs:
            valid_groups = [g for g in range(5) if group_counts[g] < 12]
            if not valid_groups:
                break
            chosen_g = random.choice(valid_groups)
            current_partition[elc] = chosen_g
            group_counts[chosen_g] += 1
            group_areas[chosen_g] += elc_areas[elc]
            
        # Verify all groups meet the 10% minimum area threshold
        if len(current_partition) == len(elc_list) and all(area >= min_area for area in group_areas.values()):
            valid_start = True

    current_score, current_metrics = evaluate_partition(
        current_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_elc
    )

    best_score = current_score
    best_metrics = current_metrics

    # Mutate and Climb
    for i in range(iterations):
        elc_to_mutate = random.choice(elc_list)
        old_group = current_partition.get(elc_to_mutate, 0)
        
        # ENFORCE MINIMUM CONSTRAINT: Cannot move if old_group will drop below 2 ER3s
        if list(current_partition.values()).count(old_group) <= 2:
            continue
            
        # ENFORCE MINIMUM AREA CONSTRAINT: Cannot move if old_group drops below 10% study area
        current_old_group_area = sum(elc_areas[e] for e, g in current_partition.items() if g == old_group)
        if (current_old_group_area - elc_areas[elc_to_mutate]) < min_area:
            continue
        
        # ENFORCE MAXIMUM CONSTRAINT: Target group must have < 12 ER3s
        valid_new_groups = []
        for g in range(5):
            if g != old_group and list(current_partition.values()).count(g) < 12:
                valid_new_groups.append(g)
                
        if not valid_new_groups:
            continue
            
        new_group = random.choice(valid_new_groups)
        
        test_partition = current_partition.copy()
        test_partition[elc_to_mutate] = new_group
        
        test_score, test_metrics = evaluate_partition(
            test_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_elc
        )
        
        if test_score >= current_score:
            current_partition = test_partition
            current_score = test_score
            
            if test_score > best_score:
                best_score = test_score
                best_metrics = test_metrics

    return best_score, best_metrics

# =====================================================================
# 2. Main Execution Block
# =====================================================================
if __name__ == '__main__':
    # Local G:\ Drive Paths
    iso_raster_path = r'G:\My Drive\Colab Notebooks\Analogs\IsoCluster.tif'
    eco3_shapefile_path = r'G:\My Drive\Colab Notebooks\Analogs\us_eco_l3_clip\us_eco_l3_clip.shp'
    zonalf_path = r'G:\My Drive\Colab Notebooks\Analogs\Zonal\Zonal_tavg.csv'

    print("Loading data and calculating master ER3 spatial overlap matrix in-memory...")
    
    raw_df = pd.read_csv(zonalf_path)
    sorted_indices = raw_df['MEAN'].sort_values(ascending=False).index
    orig_to_temp_map = {original_idx + 1: rank + 1 for rank, original_idx in enumerate(sorted_indices)}

    with rasterio.open(iso_raster_path) as src:
        iso_arr = src.read(1)
        iso_transform = src.transform

    eco_gdf = gpd.read_file(eco3_shapefile_path)
    eco_gdf['US_L3CODE'] = eco_gdf['US_L3CODE'].astype(str)
    
    # PRE-FILTER: Isolate the 15 target ER3s to prevent constraints from failing
    target_er3s = ['5', '13', '14', '18', '19', '20', '21', '22', '23', '24', '25', '26', '79', '80', '81']
    eco_gdf = eco_gdf[eco_gdf['US_L3CODE'].isin(target_er3s)]
    
    unique_elcs = eco_gdf['US_L3CODE'].unique()
    elc_to_int = {elc_str: idx + 1 for idx, elc_str in enumerate(unique_elcs)}
    int_to_elc = {idx + 1: elc_str for idx, elc_str in enumerate(unique_elcs)}
    eco_gdf['RASTER_ID'] = eco_gdf['US_L3CODE'].map(elc_to_int)

    elc_arr = rasterize(
        shapes=((geom, value) for geom, value in zip(eco_gdf.geometry, eco_gdf['RASTER_ID'])),
        out_shape=iso_arr.shape,
        transform=iso_transform,
        fill=0,
        dtype='int16'
    )

    valid_mask = (iso_arr >= 1) & (iso_arr <= 15) & (elc_arr > 0)
    iso_flat = iso_arr[valid_mask]
    elc_flat = elc_arr[valid_mask]

    overlap_matrix = pd.crosstab(iso_flat, elc_flat)
    total_valid_pixels = len(iso_flat)
    cluster_areas = pd.Series(iso_flat).value_counts()
    elc_list = list(overlap_matrix.columns)
    
    # Calculate base pixel areas for area constraint logic
    elc_areas = overlap_matrix.sum(axis=0).to_dict()
    min_area_threshold = total_valid_pixels * 0.10

    # =====================================================================
    # 3. Parallel Execution Setup
    # =====================================================================
    parallel_restarts = 48     
    iterations_per_climb = 50000 
    available_cores = multiprocessing.cpu_count()
    
    print(f"\nInitializing FORCED 100% COVERAGE ER3 Optimization...")
    print(f"Constraints: 5 Groups | Min 2/Max 12 ER3s | Min 2/Max 12 Clusters | Min 10% Area")
    print(f"Deploying {parallel_restarts} independent workers across {available_cores} CPU cores.")

    worker_args = [
        (seed, iterations_per_climb, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, elc_list, int_to_elc, elc_areas, min_area_threshold)
        for seed in range(parallel_restarts)
    ]

    global_best_score = 0
    global_best_metrics = None

    with ProcessPoolExecutor(max_workers=available_cores) as executor:
        futures = [executor.submit(run_hill_climb, arg) for arg in worker_args]
        
        for future in tqdm(as_completed(futures), total=parallel_restarts, desc="Optimizing ER3 Regions", unit="worker"):
            local_best_score, local_best_metrics = future.result()
            
            if local_best_score > global_best_score:
                global_best_score = local_best_score
                global_best_metrics = local_best_metrics

    # =====================================================================
    # 4. Output the Validated Global Maximum
    # =====================================================================
    metrics_df = pd.DataFrame(global_best_metrics)

    print("\nTable: Fully Unsupervised Global Optimum for ER3s (5 Groups | 10% Area Floor)")
    print("="*140)
    print(metrics_df.to_string(index=False))
    print(f"\nAbsolute Global Mean Case-Normalized Rate Found: {round(global_best_score, 4)}")
