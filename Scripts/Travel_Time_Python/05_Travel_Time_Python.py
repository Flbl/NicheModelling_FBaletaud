# Importation des librairies/modules
import numpy              # ... pour utiliser les classes array
import gdal               # ... pour importer/exporter raster
import datetime           # ... pour faire de compter les secondes
import funcAccessToRaster


# Debut du script
str_path_in = "D:/Florian/Documents/Universite/M1SBM/Stage/NicheModelling_FBaletaud"
str_raster_friction = str_path_in+"/eezNcGridBehrman.tif"
# list_coords_seeds = [[15443530.0,-2664001.0],[16106152.0,-2658658.0]]
list_coords_seeds = [[16060296.0,-2771872.0]]
str_path_out = "D:/Florian/Documents/Universite/M1SBM/Stage/NicheModelling_FBaletaud/Environment/Human"
str_file_out = "costNoumea"


#################################################################################
# Penser à modifier "array_friction" si frictions différentielles
# On multiplie ici par la resolution du raster car les valeurs "utiles" du
# ... raster de friction se limitent a 1. 
dict_spref = funcAccessToRaster.buildDictRasterSpatialRef(str_raster_friction)
float_res = round(dict_spref["CELLSIZE_X"],1)
raster_data = gdal.Open(str_raster_friction)
array_friction = raster_data.ReadAsArray()
array_friction = array_friction * float_res
#################################################################################


# Conversion des coordonnees en "position" pixels dans le raster a analyser
list_pixels_seeds = funcAccessToRaster.calcListCoordToPixel(list_coords_seeds, dict_spref)
                                                          

# Estimation des temps de trajets et export au format IMG
dt0 = datetime.datetime.now()
int_nb_seeds = len(list_pixels_seeds)
print "Init:", dt0
for int_id_seed_init in range(int_nb_seeds):
   print "... Seed %s" % (int_id_seed_init)
   # Calcul du raster de moindre cout
   cra = funcAccessToRaster.calcCostRasterArray(array_friction,
                                                [list_pixels_seeds[int_id_seed_init]])
   if cra == -1:
      print "... Merdouillage: funcAccessToRaster.calcCostRasterArray"
      raise error
   # Transformation des inf => nan
   cra[0][numpy.isinf(cra[0])] = -1
   # Exportation du raster de moindre cout
   int_result = funcAccessToRaster.exportCostRasterArrayToImgFiles(cra, dict_spref,
                                                                   -1, str_path_out,
                                                                   str_file_out+"_"+str(int_id_seed_init))
   if int_result == -1:
      print "... Merdouillage: funcAccessToRaster.exportCostRasterArrayToImgFiles"       
dt1 = datetime.datetime.now()
print "End:", dt1
print "Time to process:", dt1-dt0
