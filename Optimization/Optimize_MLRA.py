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
def evaluate_partition(partition_dict, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_mlra):
    metrics = []
    total_cnr_sum = 0
    
    # 6 Abstract Groups for MLRAs
    for group_id in range(6):
        group_mlra_ids = [m_id for m_id, g in partition_dict.items() if g == group_id]
        group_elc_area = overlap_matrix[group_mlra_ids].sum().sum() if group_mlra_ids else 0
        outside_total_area = total_valid_pixels - group_elc_area
        
        cluster_benefits = []
        
        # Inner Loop: Optimize Clusters (Case-Normalized Rate)
        if group_mlra_ids:
            for cluster_id in overlap_matrix.index:
                overlap_area = overlap_matrix.loc[cluster_id, group_mlra_ids].sum()
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
            # Force the top 2 clusters even if benefit is <= 0
            assigned_clusters = [c[0] for c in cluster_benefits[:2]]
        else:
            # Take the qualifying clusters, capped at 12
            assigned_clusters = qualifying_clusters[:12]
        
        cluster_area = cluster_areas[cluster_areas.index.isin(assigned_clusters)].sum() if assigned_clusters else 0
        true_positives = overlap_matrix.loc[overlap_matrix.index.isin(assigned_clusters), group_mlra_ids].sum().sum() if assigned_clusters else 0
        
        false_positives = cluster_area - true_positives
        false_negatives = group_elc_area - true_positives
        true_negatives = outside_total_area - false_positives
        
        tpr = true_positives / group_elc_area if group_elc_area > 0 else 0
        tnr = true_negatives / outside_total_area if outside_total_area > 0 else 0
        
        cnr = (tpr + tnr) / 2
        total_cnr_sum += cnr
        accuracy = (true_positives + true_negatives) / total_valid_pixels if total_valid_pixels > 0 else 0
        
        temp_orders = sorted([orig_to_temp_map[c] for c in assigned_clusters])
        original_mlra_strings = sorted([int_to_mlra[m_id] for m_id in group_mlra_ids])
        
        metrics.append({
            'Abstract Group': f"Group {group_id + 1}",
            'Optimized Clusters': str(temp_orders),
            'Optimized MLRAs': str(original_mlra_strings),
            'Prediction Accuracy': round(accuracy, 4),
            'Case-Normalized Rate': round(cnr, 4)
        })
        
    mean_cnr = total_cnr_sum / 6
    return mean_cnr, metrics

def run_hill_climb(args):
    seed, iterations, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, mlra_list, int_to_mlra = args
    random.seed(seed)

    current_partition = {}
    group_counts = {g: 0 for g in range(6)}
    available_mlras = list(mlra_list)
    random.shuffle(available_mlras)
    
    # FORCE 100% COVERAGE & MINIMUM 2 MLRAs PER GROUP
    # Step 1: Guarantee every group gets exactly 2 MLRAs to start
    for g in range(6):
        for _ in range(2):
            if available_mlras:
                mlra = available_mlras.pop()
                current_partition[mlra] = g
                group_counts[g] += 1

    # Step 2: Distribute the remainder, enforcing the Max 12 limit
    for mlra in available_mlras:
        valid_groups = [g for g in range(6) if group_counts[g] < 12]
        if not valid_groups:
            break
        chosen_g = random.choice(valid_groups)
        current_partition[mlra] = chosen_g
        group_counts[chosen_g] += 1

    current_score, current_metrics = evaluate_partition(
        current_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_mlra
    )

    best_score = current_score
    best_metrics = current_metrics

    # Mutate and Climb
    for i in range(iterations):
        mlra_to_mutate = random.choice(mlra_list)
        old_group = current_partition.get(mlra_to_mutate, 0)
        
        # ENFORCE MINIMUM CONSTRAINT: Cannot move if old_group will drop below 2 MLRAs
        if list(current_partition.values()).count(old_group) <= 2:
            continue
        
        # ENFORCE MAXIMUM CONSTRAINT: Target group must have < 12 MLRAs
        valid_new_groups = []
        for g in range(6):
            if g != old_group and list(current_partition.values()).count(g) < 12:
                valid_new_groups.append(g)
                
        if not valid_new_groups:
            continue
            
        new_group = random.choice(valid_new_groups)
        
        test_partition = current_partition.copy()
        test_partition[mlra_to_mutate] = new_group
        
        test_score, test_metrics = evaluate_partition(
            test_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, int_to_mlra
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
    mlra_shapefile_path = r'G:\My Drive\Colab Notebooks\Analogs\mlra_clip\mlra_clip.shp'
    zonalf_path = r'G:\My Drive\Colab Notebooks\Analogs\Zonal\Zonal_tavg.csv'

    print("Loading data and calculating master MLRA spatial overlap matrix in-memory...")
    
    raw_df = pd.read_csv(zonalf_path)
    sorted_indices = raw_df['MEAN'].sort_values(ascending=False).index
    orig_to_temp_map = {original_idx + 1: rank + 1 for rank, original_idx in enumerate(sorted_indices)}

    with rasterio.open(iso_raster_path) as src:
        iso_arr = src.read(1)
        iso_transform = src.transform

    mlra_gdf = gpd.read_file(mlra_shapefile_path)
    
    # PRE-FILTER: Isolate the 33 target MLRAs to prevent constraints from failing
    target_mlras = ['30', '40', '24', '27', '28A', '29', '34A', '34B', '35', '42B',
                    '11', '13', '22A', '23', '25', '26', '28B', '36', '46', '47', '48A', '51',
                    '38', '41', '39', '42A', '42C', '70A', '70B', '77B', '77C', '77D', '77E']
    
    mlra_gdf = mlra_gdf[mlra_gdf['MLRARSYM'].isin(target_mlras)]
    
    unique_mlras = mlra_gdf['MLRARSYM'].unique()
    mlra_to_int = {mlra_str: idx + 1 for idx, mlra_str in enumerate(unique_mlras)}
    int_to_mlra = {idx + 1: mlra_str for idx, mlra_str in enumerate(unique_mlras)}
    mlra_gdf['RASTER_ID'] = mlra_gdf['MLRARSYM'].map(mlra_to_int)

    mlra_arr = rasterize(
        shapes=((geom, value) for geom, value in zip(mlra_gdf.geometry, mlra_gdf['RASTER_ID'])),
        out_shape=iso_arr.shape,
        transform=iso_transform,
        fill=0,
        dtype='int16'
    )

    valid_mask = (iso_arr >= 1) & (iso_arr <= 15) & (mlra_arr > 0)
    iso_flat = iso_arr[valid_mask]
    mlra_flat = mlra_arr[valid_mask]

    overlap_matrix = pd.crosstab(iso_flat, mlra_flat)
    total_valid_pixels = len(iso_flat)
    cluster_areas = pd.Series(iso_flat).value_counts()
    mlra_list = list(overlap_matrix.columns)

    # =====================================================================
    # 3. Parallel Execution Setup
    # =====================================================================
    parallel_restarts = 50     
    iterations_per_climb = 50000 
    available_cores = multiprocessing.cpu_count()
    
    print(f"\nInitializing FORCED 100% COVERAGE MLRA Optimization (Min 2/Max 12 MLRAs | Min 2/Max 12 Clusters)...")
    print(f"Deploying {parallel_restarts} independent workers across {available_cores} CPU cores.")

    worker_args = [
        (seed, iterations_per_climb, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, mlra_list, int_to_mlra)
        for seed in range(parallel_restarts)
    ]

    global_best_score = 0
    global_best_metrics = None

    # Implement tqdm progress bar with as_completed
    with ProcessPoolExecutor(max_workers=available_cores) as executor:
        futures = [executor.submit(run_hill_climb, arg) for arg in worker_args]
        
        for future in tqdm(as_completed(futures), total=parallel_restarts, desc="Optimizing MLRA Regions", unit="worker"):
            local_best_score, local_best_metrics = future.result()
            
            if local_best_score > global_best_score:
                global_best_score = local_best_score
                global_best_metrics = local_best_metrics

    # =====================================================================
    # 4. Output the Validated Global Maximum
    # =====================================================================
    metrics_df = pd.DataFrame(global_best_metrics)

    print("\nTable: Fully Unsupervised Global Optimum for MLRAs (Min 2 & Max 12 Constraints Forced)")
    print("="*140)
    print(metrics_df.to_string(index=False))
    print(f"\nAbsolute Global Mean Case-Normalized Rate Found: {round(global_best_score, 4)}")
    
