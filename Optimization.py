import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
import random
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import warnings
warnings.filterwarnings('ignore')

# =====================================================================
# 1. Top-Level Worker Functions 
# =====================================================================
def evaluate_partition(partition_dict, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map):
    metrics = []
    total_cnr_sum = 0
    
    for group_id in range(5):
        # Only evaluate ELCs actively assigned to groups 0-4 (Ignore Unassigned pool -1)
        group_er3s = [er3 for er3, g in partition_dict.items() if g == group_id]
        group_elc_area = overlap_matrix[group_er3s].sum().sum() if group_er3s else 0
        outside_total_area = total_valid_pixels - group_elc_area
        
        qualifying_clusters = []
        
        # Inner Loop: Optimize Clusters (Case-Normalized Rate)
        if group_er3s:
            for cluster_id in overlap_matrix.index:
                overlap_area = overlap_matrix.loc[cluster_id, group_er3s].sum()
                outside_spill = cluster_areas.get(cluster_id, 0) - overlap_area
                
                tpr_gain = overlap_area / group_elc_area if group_elc_area > 0 else 0
                tnr_loss = outside_spill / outside_total_area if outside_total_area > 0 else 0
                
                benefit = tpr_gain - tnr_loss
                if benefit > 0:
                    qualifying_clusters.append((cluster_id, benefit))
                    
        # ENFORCE CONSTRAINT: Max 7 clusters per group (Sort by highest benefit)
        qualifying_clusters.sort(key=lambda x: x[1], reverse=True)
        assigned_clusters = [c[0] for c in qualifying_clusters[:7]]
        
        # Calculate Group Metrics
        cluster_area = cluster_areas[cluster_areas.index.isin(assigned_clusters)].sum() if assigned_clusters else 0
        true_positives = overlap_matrix.loc[overlap_matrix.index.isin(assigned_clusters), group_er3s].sum().sum() if assigned_clusters else 0
        
        false_positives = cluster_area - true_positives
        false_negatives = group_elc_area - true_positives
        true_negatives = outside_total_area - false_positives
        
        tpr = true_positives / group_elc_area if group_elc_area > 0 else 0
        tnr = true_negatives / outside_total_area if outside_total_area > 0 else 0
        
        cnr = (tpr + tnr) / 2
        total_cnr_sum += cnr
        accuracy = (true_positives + true_negatives) / total_valid_pixels if total_valid_pixels > 0 else 0
        
        temp_orders = sorted([orig_to_temp_map[c] for c in assigned_clusters])
        metrics.append({
            'Abstract Group': f"Group {group_id + 1}",
            'Optimized Clusters': str(temp_orders),
            'Optimized ER3s': str(sorted(group_er3s)),
            'Prediction Accuracy': round(accuracy, 4),
            'Case-Normalized Rate': round(cnr, 4)
        })
        
    mean_cnr = total_cnr_sum / 5
    return mean_cnr, metrics

def run_hill_climb(args):
    seed, iterations, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, er3_list = args
    random.seed(seed)

    # Initialize with -1 (Unassigned Pool) to prevent accidental super-groups
    current_partition = {er3: -1 for er3 in er3_list}
    available_er3s = list(er3_list)
    random.shuffle(available_er3s)
    
    for g in range(5):
        num_to_assign = random.randint(1, 7)
        for _ in range(num_to_assign):
            if available_er3s:
                current_partition[available_er3s.pop()] = g

    current_score, current_metrics = evaluate_partition(
        current_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map
    )

    best_score = current_score
    best_metrics = current_metrics

    # Mutate and Climb
    for i in range(iterations):
        er3_to_mutate = random.choice(er3_list)
        old_group = current_partition[er3_to_mutate]
        
        # Determine valid moves (Group -1 is always valid, Groups 0-4 are valid if they have < 7 ELCs)
        valid_new_groups = [-1]
        for g in range(5):
            if g != old_group and list(current_partition.values()).count(g) < 7:
                valid_new_groups.append(g)
                
        if not valid_new_groups:
            continue
            
        new_group = random.choice(valid_new_groups)
        
        test_partition = current_partition.copy()
        test_partition[er3_to_mutate] = new_group
        
        test_score, test_metrics = evaluate_partition(
            test_partition, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map
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

    print("Loading data and calculating master spatial overlap matrix in-memory...")
    
    raw_df = pd.read_csv(zonalf_path)
    sorted_indices = raw_df['MEAN'].sort_values(ascending=False).index
    orig_to_temp_map = {original_idx + 1: rank + 1 for rank, original_idx in enumerate(sorted_indices)}

    with rasterio.open(iso_raster_path) as src:
        iso_arr = src.read(1)
        iso_transform = src.transform

    elc_gdf = gpd.read_file(eco3_shapefile_path)
    elc_gdf['US_L3CODE'] = elc_gdf['US_L3CODE'].astype(int)

    elc_arr = rasterize(
        shapes=((geom, value) for geom, value in zip(elc_gdf.geometry, elc_gdf['US_L3CODE'])),
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
    er3_list = list(overlap_matrix.columns)

    # =====================================================================
    # 3. Parallel Execution Setup
    # =====================================================================
    parallel_restarts = 30     
    iterations_per_climb = 25000 
    available_cores = multiprocessing.cpu_count()
    
    print(f"\nInitializing Constrained Optimization (Max 7 ELCs & 7 Clusters per group)...")
    print(f"Deploying {parallel_restarts} independent workers across {available_cores} CPU cores.")

    worker_args = [
        (seed, iterations_per_climb, overlap_matrix, total_valid_pixels, cluster_areas, orig_to_temp_map, er3_list)
        for seed in range(parallel_restarts)
    ]

    global_best_score = 0
    global_best_metrics = None

    with ProcessPoolExecutor(max_workers=available_cores) as executor:
        results = executor.map(run_hill_climb, worker_args)
        
        for local_best_score, local_best_metrics in results:
            if local_best_score > global_best_score:
                global_best_score = local_best_score
                global_best_metrics = local_best_metrics

    # =====================================================================
    # 4. Output the Validated Global Maximum
    # =====================================================================
    metrics_df = pd.DataFrame(global_best_metrics)

    print("\nTable: Fully Unsupervised Global Optimum (Constrained: Max 7 ELCs/Clusters)")
    print("="*125)
    print(metrics_df.to_string(index=False))
    print(f"\nAbsolute Global Mean Case-Normalized Rate Found: {round(global_best_score, 4)}")
