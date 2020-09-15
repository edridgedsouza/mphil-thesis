library(ClusterMap)
setwd('/home/edsouza/Documents/Research/RussellLab/')

marker_file_list <- c(zinzen = '2020_04_29_ClusterMap/zinzen.markers.csv', 
                      avalos = '2020_04_29_ClusterMap/avalos.markers.csv')

single_obj_list <- c(zinzen = zinzen, avalos = avalos)

# res <- cluster_map(marker_file_list, edge_cutoff = 0.1, 
#                    output = '2020_04_29_ClusterMap/ZinzenAvalos', 
#                    single_obj_list = single_obj_list)


cell_num_list <- lapply(single_obj_list,
                        function(obj){
                          summary(Idents(obj))
                        })
mapRes <- cluster_map_by_marker(marker_file_list, cutoff = 0.1, output = '2020_04_29_ClusterMap/ZinzenAvalos')
circos_map(mapRes, cell_num_list,  output = '2020_04_29_ClusterMap/ZinzenAvalos')

# Found this info by hacking through the debug before the circos plot was made
# Browse[3]> pair
# v1        v2 similarity regroup
# 1   zinzen_5  avalos_4       0.09       1
# 2   zinzen_5  avalos_5       0.09       1
# 3  zinzen_10 avalos_10       0.04       3
# 4  zinzen_10 avalos_11       0.04       3
# 5  zinzen_10 avalos_12       0.04       3
# 6  zinzen_10 avalos_13       0.04       3
# 7  zinzen_10 avalos_14       0.04       3
# 8  zinzen_10 avalos_15       0.04       3
# 9  zinzen_10 avalos_16       0.04       3
# 10 zinzen_10 avalos_17       0.04       3
# 11 zinzen_10 avalos_18       0.04       3
# 12 zinzen_10  avalos_6       0.04       3
# 13 zinzen_10  avalos_7       0.04       3
# 14 zinzen_10  avalos_8       0.04       3