# !/usr/bin/env python
# -*- coding: utf-8 -*


import datetime
import gdal
import numpy
import os
import psutil
import sys
from skimage import graph


NUMPY_TO_GDAL = {numpy.dtype("int8"):None,
                 numpy.dtype("int16"):gdal.GDT_Int16,
                 numpy.dtype("int32"):gdal.GDT_Int32,
                 numpy.dtype("int64"):None,
                 numpy.dtype("uint8"):gdal.GDT_Byte,
                 numpy.dtype("uint16"):gdal.GDT_UInt16,
                 numpy.dtype("uint32"):gdal.GDT_UInt32,
                 numpy.dtype("uint64"):None,
                 numpy.dtype("float16"):None,
                 numpy.dtype("float32"):gdal.GDT_Float32,
                 numpy.dtype("float64"):gdal.GDT_Float64}
                  

def buildDictIndexFirstAndLastValues(array, var_invalid_value):
    """
        Construit un dictionnaire stockant les indices de la premiere et derniere
        valeur valide
    """
    try:
        dict_idx_fl  = {"IDX_INI":None, "IDX_END":None}
        int_nb_values = len(array)
        array_bool = numpy.where(array == var_invalid_value, False, True)
        int_idx_ini = numpy.argmax(numpy.cumsum(array_bool) > 0)
        dict_idx_fl["IDX_INI"] = int_idx_ini
        int_idx_end = int_nb_values - numpy.argmax(numpy.cumsum(array_bool[::-1]) > 0)
        dict_idx_fl["IDX_END"] = int_idx_end - 1
        return dict_idx_fl
    except:
        return -1
    
def buildDictRasterSpatialRef(str_raster):
    """
        Construit un dictionnaire stockant la reference spatiale d'un Raster.
    """
    try:
        dict_spref = {"CELLSIZE_X":None,"CELLSIZE_Y":None,
                      "ULCORNER_X":None,"ULCORNER_Y":None,
                      "ULPIXEL_COL":None,"ULPIXEL_ROW":None,
                      "N_COLS":None,"N_ROWS":None}
        raster_in = gdal.Open(str_raster)
        tupple_info = raster_in.GetGeoTransform()
        int_nb_cols = int(raster_in.RasterXSize)
        int_nb_rows = int(raster_in.RasterYSize)
        dict_spref["CELLSIZE_X"] = abs(tupple_info[1])
        dict_spref["CELLSIZE_Y"] = abs(tupple_info[5])
        dict_spref["ULCORNER_X"] = tupple_info[0]
        dict_spref["ULCORNER_Y"] = tupple_info[3]
        dict_spref["ULPIXEL_COL"] = tupple_info[2]
        dict_spref["ULPIXEL_ROW"] = tupple_info[4]
        dict_spref["N_COLS"] = int_nb_cols
        dict_spref["N_ROWS"] = int_nb_rows
        return dict_spref
    except:
        return -1

def calcListCoordToPixel(list_coords, dict_spref):
    """
        Calcule le positionnement image d'une liste de coordonnees a partir
        du dictionnaire de reference spatiale d'un Raster. 
    """
    try:
        int_nb_coords = len(list_coords)
        array_coords = numpy.array(list_coords, dtype=numpy.float64)
        array_coords = numpy.transpose(array_coords, (1,0))
        float_origin_x = dict_spref["ULCORNER_X"]
        if dict_spref["ULPIXEL_COL"] != 0:
            float_origin_x -= dict_spref["ULPIXEL_COL"]*dict_spref["CELLSIZE_X"]
        float_origin_y = dict_spref["ULCORNER_Y"]
        if dict_spref["ULPIXEL_ROW"] != 0:
            float_origin_y += dict_spref["ULPIXEL_ROW"]*dict_spref["CELLSIZE_Y"]
        array_coords[0] = (array_coords[0] - float_origin_x) / dict_spref["CELLSIZE_X"]
        array_coords[0] = numpy.trunc(array_coords[0])
        array_coords[1] = (float_origin_y - array_coords[1]) / dict_spref["CELLSIZE_Y"]
        array_coords[1] = numpy.trunc(array_coords[1])
        array_coords = numpy.transpose(array_coords, (1,0))
        array_coords = array_coords.astype(int)
        list_pixels = [list(array_coords[i]) for i in range(array_coords.shape[0])]
        return list_pixels
    except:
        raise error


def calcListPixelToCoord(list_pixels, dict_spref, str_type="Center"):
    """
        Calcule les coordonnees d'une liste de position pixel a partir
        du dictionnaire de reference spatiale d'un Raster. 
    """
    try:
        int_nb_pixels = len(list_pixels)
        array_pixels = numpy.array(list_pixels, dtype=numpy.float64)
        array_pixels = numpy.transpose(array_pixels, (1,0))
        if str_type == "Center":
            float_add = 0.5
        elif str_type == "UpperLeft":
            float_add = 0.0
        else:
            raise error
        float_origin_x = dict_spref["ULCORNER_X"]
        if dict_spref["ULPIXEL_COL"] != 0:
            float_origin_x -= dict_spref["ULPIXEL_COL"]*dict_spref["CELLSIZE_X"]
        float_origin_y = dict_spref["ULCORNER_Y"]
        if dict_spref["ULPIXEL_ROW"] != 0:
            float_origin_y += dict_spref["ULPIXEL_ROW"]*dict_spref["CELLSIZE_Y"]
        array_pixels[0] = float_origin_x + (array_pixels[0] + float_add) * dict_spref["CELLSIZE_X"]
        array_pixels[1] = float_origin_y - (array_pixels[1] + float_add) * dict_spref["CELLSIZE_Y"]
        array_pixels = numpy.transpose(array_pixels, (1,0))
        list_coords = [list(array_pixels[i]) for i in range(array_pixels.shape[0])]
        return list_coords
    except:
        raise error

def calcListPixelWithNewOrigin(list_pixels, list_pixel_origin):
    """
        Calcule le positionnement d'une liste de pixel a partir d'une nouvelle origine. 
    """
    try:
        int_nb_pixels = len(list_pixels)
        array_pixels = numpy.array(list_pixels)
        array_pixel_origin = numpy.array(list_pixel_origin)
        array_pixels = array_pixels - array_pixel_origin
        list_pixels = [list(array_pixels[i]) for i in range(array_pixels.shape[0])]
        return list_pixels
    except:
        raise error

def calcListExtentForAccessTo(list_coords, dict_spref, float_access_max, float_friction_min):
    """
        Calcule l'etendue necessaire pour une analyse AccessTo.
    """
    try:
        list_extent = [None, None, None, None]
        int_nb_coord = len(list_coords)
        float_radius = (float(float_access_max) / float_friction_min) + 1
        float_radius_x = float_radius * dict_spref["CELLSIZE_X"]
        float_radius_y = float_radius * dict_spref["CELLSIZE_Y"]
        list_eflc = calcListExtentFromListCoords(list_coords)
        if list_eflc == -1:
            raise error
        if int_nb_coord == 1: 
            list_extent = [[list_eflc[0][0]-float_radius_x, list_eflc[0][1]+float_radius_y],
                           [list_eflc[0][0]+float_radius_x, list_eflc[0][1]-float_radius_y]]     
        else:
            list_extent = [[list_eflc[0][0]-float_radius_x, list_eflc[0][1]+float_radius_y],
                           [list_eflc[1][0]+float_radius_x, list_eflc[1][1]-float_radius_y]]
        return list_extent
    except:
        return -1

def calcListExtentFromListCoords(list_coords):
    """
        Calcule l'etendue d'une liste des coordonnees.
    """
    try:
        array_coords = numpy.array(list_coords, dtype=numpy.float64)
        array_coords_min = numpy.min(array_coords, 0)
        array_coords_max = numpy.max(array_coords, 0)
        list_extent = [[array_coords_min[0], array_coords_max[1]],
                       [array_coords_max[0], array_coords_min[1]]]
        return list_extent
    except:
        return -1

################################################################################
######   VOIR COMMENT AMELIORER / COMPLETER LA FONCTION AVEC DES GRAPHS   ######
######                     AVEC PROPRIETES DIFFERENTES                    ######
def calcCostRasterArray(array, list_pixels):
    """
        Calcul le raster de moindre cout geometrique.
    """
    try:
        graph_mcp = graph.MCP_Geometric(array)
        list_pixels_rev = [[list_pixels[i][1],list_pixels[i][0]]
                           for i in range(len(list_pixels))]
        array_mcp = graph_mcp.find_costs(list_pixels_rev)
        array_mcp = [array_mcp[0], array_mcp[1]]
        return array_mcp
    except:
        return -1
################################################################################
################################################################################
################################################################################


def exportCostRasterArrayToImgFiles(array_mcp, dict_spref, var_no_data, str_path,
                                    str_file_cost="", str_file_link=""):
    """
        Exporte les array de moindre cout geometrique sur le disque au format imagine.
    """
    try:
        tupple_size = array_mcp[0].shape
        int_nb_cols = tupple_size[1]
        int_nb_rows = tupple_size[0]
        list_geo_transform = [dict_spref["ULCORNER_X"], dict_spref["CELLSIZE_X"],
                              dict_spref["ULPIXEL_COL"], dict_spref["ULCORNER_Y"],
                              dict_spref["ULPIXEL_ROW"], -dict_spref["CELLSIZE_Y"]]
        driver = gdal.GetDriverByName("HFA")
        if str_file_cost != "":
            str_path_file_cost = str_path+"/"+str_file_cost+".img"
            type_gdal = NUMPY_TO_GDAL[array_mcp[0].dtype]
            raster_out = driver.Create(str_path_file_cost, int_nb_cols, int_nb_rows, 1, type_gdal)
            raster_out.SetGeoTransform(list_geo_transform)
            raster_out_band = raster_out.GetRasterBand(1)
            raster_out_band.SetNoDataValue(var_no_data)
            raster_out_band.WriteArray(array_mcp[0])
            del raster_out
        return 1
    except:
        return -1

#def expandRasterArray(array, ...):
#    """
#        Agrandit l'etendue d'un array.
#    """
#    try:
#        print 45
#    except:
#        return -1



def shrinkRasterArray(array, str_type, list_params, dict_spref={}):
    """
        Reduit l'etendue d'un array.
    """
    try:
        dict_extent = {"ROW_INI":None,"ROW_END":None,
                       "COL_INI":None,"COL_END":None}
        if str_type == "ByValue":
            array_mask = numpy.ones(array.shape, numpy.uint8)
            for var_param in list_params:
                array_mask = numpy.where(array == var_param, 0, array_mask)
            dict_indices = {0:None, 1:None}
            for int_id_axe in [0,1]:
                array_sum = numpy.sum(array_mask, int_id_axe)
                dict_ifalv = buildDictIndexFirstAndLastValues(array_sum, 0)
                if dict_ifalv == -1:
                    raise error
                dict_indices[int_id_axe] = dict_ifalv
            dict_extent["COL_INI"] = dict_indices[0]["IDX_INI"]
            dict_extent["COL_END"] = dict_indices[0]["IDX_END"] + 1
            dict_extent["ROW_INI"] = dict_indices[1]["IDX_INI"]
            dict_extent["ROW_END"] = dict_indices[1]["IDX_END"] + 1
        elif str_type == "ByExtent":
            dict_extent["COL_INI"] = list_params[0][0]
            dict_extent["COL_END"] = list_params[1][0] + 1
            dict_extent["ROW_INI"] = list_params[0][1]
            dict_extent["ROW_END"] = list_params[1][1] + 1
        else:
            raise error
        if dict_spref != {}:
            dict_spref = shrinkDictRasterSpatialRef(dict_spref, dict_extent)
            if dict_spref == -1:
                raise error
        dict_results = {"ROW_INI":dict_extent["ROW_INI"],
                        "COL_INI":dict_extent["COL_INI"],
                        "DICT_SPREF":dict_spref,
                        "ARRAY":array[dict_extent["ROW_INI"]:dict_extent["ROW_END"],
                                      dict_extent["COL_INI"]:dict_extent["COL_END"]]}
        return dict_results
    except:
        return -1

def shrinkDictRasterSpatialRef(dict_spref, dict_extent):
    """
        Met a jour le dictionnaire stockant la reference spatiale d'un Raster reduit.
    """
    try:
        dict_spref_shrink = dict_spref.copy()
        dict_spref_shrink["N_COLS"] = int(dict_extent["COL_END"] - dict_extent["COL_INI"])
        dict_spref_shrink["N_ROWS"] = int(dict_extent["ROW_END"] - dict_extent["ROW_INI"])
        dict_spref_shrink["ULPIXEL_COL"] = 0
        dict_spref_shrink["ULPIXEL_ROW"] = 0
        int_delta_col = dict_extent["COL_INI"] - dict_spref["ULPIXEL_COL"]
        float_coord_x_origin = dict_spref["ULCORNER_X"] + dict_spref["CELLSIZE_X"] * int_delta_col
        dict_spref_shrink["ULCORNER_X"] = float_coord_x_origin
        int_delta_row = dict_extent["ROW_INI"] - dict_spref["ULPIXEL_ROW"]
        float_coord_y_origin = dict_spref["ULCORNER_Y"] - dict_spref["CELLSIZE_Y"] * int_delta_row
        dict_spref_shrink["ULCORNER_Y"] = float_coord_y_origin
        return dict_spref_shrink
    except:
        return -1

def testListPathsNok(list_paths):
    """
        Test si plusieurs chemins sont absent du disque.  
    """
    try:
        list_paths_nok = []
        for str_path in list_paths:
            if os.path.exists(str_path) == False:
                list_paths_nok.append(str_path)
        return list_paths_nok
    except:
        return -1

def testListRasterSpatialRef(list_rasters):
    """
        Test si une plusieurs Raster ont la meme reference spatiale.
    """
    try:
        dict_spref = {"CELLSIZE_X":None,"CELLSIZE_Y":None,
                      "ULCORNER_X":None,"ULCORNER_Y":None,
                      "ULPIXEL_COL":None,"ULPIXEL_ROW":None,
                      "N_COLS":None,"N_ROWS":None}
        list_keys_spref = dict_spref.keys()
        for list_raster in list_rasters:
            dict_rsr = buildDictRasterSpatialRef(list_raster[0])
            if dict_rsr == -1:
                raise error
            for key in list_keys_spref:
                if dict_spref[key] == None:
                    dict_spref[key] = dict_rsr[key]
                elif dict_spref[key] != dict_rsr[key]:
                    dict_spref[key] = False
        return dict_spref
    except:
        return -1
