
import arcpy, os
import arcpy.sa

# Set workspace, scratch workspace, and output here
# Note: Workspace and Scratch can be same thing if needed
workspace = r"C:\Users\chris\Documents\ArcGIS\Projects\MEPModelCombo\Input.gdb"
scratch = r"C:\Users\chris\Documents\ArcGIS\Projects\MEPModelCombo\Intermediate.gdb"
output = r"C:\Users\chris\Documents\ArcGIS\Projects\MEPModelCombo\Output.gdb"

# Do you want intermediate layers to be deleted? If so make variable true.
delete_extras_wanted = True

# Layer Names
polylines = "landDx_polylines"
clip = "rsf_region_zones"
albedo = "Albedo_2011_2017_from2001_2010"
points = "landDx_points"

### Layers for LULC
bare_built_rock = "bare_built_rock"
crop = "crop"
degraded = "degraded"
agriculture = "Mara_Agriculture"

# Wanted Projection
projection = "PROJCS['WGS_1984_UTM_Zone_36S',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984'," \
             "SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0]," \
             "UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator']," \
             "PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',10000000.0]," \
             "PARAMETER['Central_Meridian',33.0],PARAMETER['Scale_Factor',0.9996]," \
             "PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"

# Cell Size
cellsize = 30

########################################################################################################################
'''***Change the weights or classifications of the attributes of each layer below***'''
# *** FENCE TYPE WEIGHTS ***
wire_fence_type = 0.3
electric_fence_type = 0.5
other_fence_type = 0.2

# *** LAND OWNERSHIP TYPE CLASSIFICATION ***
protected_land_type = 1
community_conservancy_land_type = 3
unprotected_land_type = 5

# *** LAND USE/LAND CHANGE CLASSIFICATION ***
bare_rock_LCLU = 1
crop_LCLU = 5
degraded_land_LCLU = 2
agriculture_LCLU = 5


# *** ROAD DISTANCE CLASSIFICATION ****
# In each bracket the first 2 numbers are the range of distances and the last value is the value
# EX: [0, 250, 5] means anything within the distance of 0 - 250 meters to road type is going to have a value of 5
primary_road_type = [[0, 250, 5], [250, 500, 4], [500, 750, 3], [750, 1000, 2], [1000, 5000, 1], [5000, 10000000, 0]]
secondary_road_type = [[0, 250, 3], [250, 500, 2], [500, 750, 1], [750, 1000, 1], [1000, 5000, 1], [5000, 10000000, 0]]
tertiary_road_type = [[0, 250, 2], [250, 500, 1], [500, 750, 1], [750, 1000, 0], [1000, 5000, 0], [5000, 10000000, 0]]

# *** SETTLEMENT DISTANCE CLASSIFICATION ****
# In each bracket the first 2 numbers are the range of distances and the last value is the value
# EX: [0, 400, 1] means anything within the distance of 0 - 400 meters to settlement type is going to have a value of 1
boma_settlement_type = [[0, 400, 1], [400, 800, 1], [800, 1200, 0], [1200, 1600, 0], [1600, 5000, 0], [5000, 10000000, 0]]
town_settlement_type = [[0, 400, 5], [400, 800, 4], [800, 1200, 3], [1200, 1600, 2], [1600, 5000, 1], [5000, 10000000, 0]]
safari_settlement_type = [[0, 400, 3], [400, 800, 3], [800, 1200, 2], [1200, 1600, 1], [1600, 5000, 0], [5000, 10000000, 0]]
hotel_settlement_type = [[0, 400, 4], [400, 800, 4], [800, 1200, 3], [1200, 1600, 2], [1600, 5000, 1], [5000, 10000000, 0]]

# **** FINAL LAYER WEIGHTS ***
# *** Change final weighted sum weights for all layers here! ***
Fence_Weight = .10
LCLU_Weight = .20
Land_Ownership_Weight = .10
Roads_Weight = .20
Albedo_Weight = .05
Population_Density_Weight = .10
Settlement_Weight = .25

#######################################################################################################################
''' ***Below are the functions for each of the wanted layers. Instructions for changing the classifications and 
weights are above***'''

def FenceModel(wire_fence_weight, electric_fence_weight, other_fence_weight):
    print("Starting Fence Processing:")
    delete_list = []
    arcpy.env.overwriteOutput = True
    # Process: Project (Project) (management)
    with arcpy.EnvManager(scratchWorkspace = scratch, workspace = workspace):
        print("     Projecting and Clipping")
        # Process: Project (Project) (management)
        out_data = polylines + "UTM"
        out_data = os.path.join(scratch, out_data)
        arcpy.management.Project(polylines, out_data, projection)
        delete_list.append(out_data)
        # Process: Clip (Clip) (analysis)
        out_data1 = polylines + "UTM_Clip"
        out_data1 = os.path.join(scratch, out_data1)
        arcpy.analysis.Clip(out_data, clip, out_data1)
        delete_list.append(out_data1)
        extent = arcpy.Describe(clip).extent
        print("     Selecting layers by fence type.")
        # Process: Make Electric Fence Layer
        electric_fence = "Electric_Fence"
        where_clause = "type LIKE '%Fence\\_Electric%' ESCAPE '\\'"
        arcpy.management.MakeFeatureLayer(out_data1, electric_fence, where_clause)

        # Process: Make Wire Fence Layer
        wire_fence = "Wire_Fence"
        where_clause="type LIKE '%Fence\\_Wire%' ESCAPE '\\'"
        arcpy.management.MakeFeatureLayer(out_data1, wire_fence, where_clause)
        # Process: Make Other Fence Layer
        other_fence = "Other_Fence"
        where_clause = """type NOT LIKE '%Fence\\_Electric%' ESCAPE '\\' And type NOT LIKE '%Fence\\_Wire%' ESCAPE '\\' And
        type LIKE '%Fence%' ESCAPE '\\'"""
        arcpy.management.MakeFeatureLayer(out_data1, other_fence, where_clause)
        print("     Creating line density layers.")
        # Electric Fence Line Density
        with arcpy.EnvManager(extent=extent):
            EF_density = os.path.join(scratch, "electric_fence_density")
            EF_density_tool = arcpy.sa.LineDensity(electric_fence, population_field="NONE",
                                                      cell_size=cellsize, search_radius=1609.334,
                                                      area_unit_scale_factor="SQUARE_METERS")
            EF_density_tool.save(EF_density)
            delete_list.append(EF_density)
        # Wire Fence Line Density
        with arcpy.EnvManager(extent=extent):
            WF_density = os.path.join(scratch, "wire_fence_density")
            WF_density_tool = arcpy.sa.LineDensity(wire_fence, population_field="NONE",
                                          cell_size=cellsize, search_radius=1609.334,
                                          area_unit_scale_factor="SQUARE_METERS")
            WF_density_tool.save(WF_density)
            delete_list.append(WF_density)
        # Other Fence Line Density
        with arcpy.EnvManager(extent=extent):
            OF_density = os.path.join(scratch, "other_fence_density")
            OF_density_tool = arcpy.sa.LineDensity(other_fence, population_field="NONE",
                                          cell_size=cellsize, search_radius=1609.334,
                                          area_unit_scale_factor="SQUARE_METERS")
            OF_density_tool.save(OF_density)
            delete_list.append(OF_density)
        print("     Summing line density layers.")
        # Weighted Sum
        weighted_fence = os.path.join(scratch, "weighted_fence")
        with arcpy.EnvManager(extent="MINOF", mask=clip, scratchWorkspace= scratch,workspace= workspace):

            WSumTableObj = arcpy.sa.WSTable([[WF_density, "VALUE", wire_fence_weight],
                                            [OF_density, "VALUE", other_fence_weight],
                                             [EF_density, "VALUE", electric_fence_weight]])
            weighted_fence_tool = arcpy.sa.WeightedSum(WSumTableObj)
            weighted_fence_tool.save(weighted_fence)
            delete_list.append(weighted_fence)

        print("     Rescaling to 0-5.")
        # Rescale by function
        scaled_fenceWeight_FINAL = os.path.join(scratch, "FenceWeight_almostFINAL")
        min = arcpy.sa.Raster(weighted_fence).minimum
        max = arcpy.sa.Raster(weighted_fence).maximum
        scaled_fenceweight_tool = arcpy.sa.RescaleByFunction(in_raster= weighted_fence,
                                                             transformation_function=arcpy.sa.TfLogisticDecay(min, max),
                                                             from_scale=5, to_scale=0.0)
        scaled_fenceweight_tool.save(scaled_fenceWeight_FINAL)

        delete_list.append(scaled_fenceWeight_FINAL)
        with arcpy.EnvManager(mask=clip):
            out_raster = arcpy.sa.Reclassify(scaled_fenceWeight_FINAL,
                                             "VALUE", "0 0;0 1 1;1 2 2;2 3 3;3 4 4;4 5 5;NODATA 0", "DATA")
        finalName = "FenceWeight_FINAL"
        fenceWeight_FINAL = os.path.join(output, finalName)
        out_raster.save(fenceWeight_FINAL)
        FenceModel.name = fenceWeight_FINAL

    print(f"Fence Layer Created!\n")
    delete_extras(delete_list)
    return finalName


def Land_Ownership(mr_weight, cca_weight, up_weight):
    print("Starting Land Ownership Processing:")
    ownership_layer = clip
    delete_list = []
    arcpy.env.overwriteOutput = True
    with arcpy.EnvManager(scratchWorkspace=scratch, workspace=workspace):
        print("     Projecting")
    # Process: Project (Project) (management)
        out_data = "GME_Extent_UTM_36S"
        out_data = os.path.join(scratch, out_data)
        arcpy.management.Project(ownership_layer, out_data, projection)
        delete_list.append(out_data)
        print("     Converting Polygon to Raster")
    # Process: Polygon to Raster (Polygon to Raster) (conversion)
        Poly_to_Raster = out_data + "_raster"
        arcpy.conversion.PolygonToRaster(in_features=out_data, value_field="zone",
                                         out_rasterdataset=Poly_to_Raster, cell_assignment="CELL_CENTER",
                                         priority_field="NONE", cellsize=cellsize, build_rat="BUILD")
        print("     Reclassifying")
        delete_list.append(Poly_to_Raster)
    # Process: Reclassify (Reclassify) (sa)
    finalName = "LandOwnership_FINAL"
    outRaster = os.path.join(output, finalName)

    ownership_remap = f"mr {mr_weight};cca {cca_weight};up {up_weight}"
    GME_ownership_36S_Reclassify = arcpy.sa.Reclassify(in_raster=Poly_to_Raster, reclass_field="zone",
                                                       remap=ownership_remap, missing_values="NODATA")
    GME_ownership_36S_Reclassify.save(outRaster)
    Land_Ownership.name = outRaster
    delete_extras(delete_list)
    print(f"Land Ownership Layer Complete!\n")
    return finalName

def LULC(bare_rock_weight, crop_weight, degraded_weight, agriculture_weight):
    print(f"Starting LULC processing: ")
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = workspace
    rasterlist = [bare_built_rock, crop, degraded]
    newlist = []
    print("     Projecting and Clipping rasters.")
    for x in rasterlist:
        name = x+"_projMask"
        path = os.path.join(scratch, name)
        arcpy.management.ProjectRaster(arcpy.Raster(x), path, projection)
        outExtractMask = arcpy.sa.ExtractByMask(path, clip)
        outExtractMask.save(path)
        newlist.append(path)
    name = agriculture+"_proj"
    path = os.path.join(scratch, name)
    maskpath = os.path.join(scratch, name+"Mask")
    arcpy.management.Project(agriculture, path, projection)
    arcpy.analysis.Clip(path, clip, maskpath)
    arcpy.management.Delete(path)
# Reclassify
    print(f"     Reclassifying")
    delete_list = []
    for x in newlist:
        if "bare_built_rock" in x:
            Reclass_bare1 = x
            Reclassify_3_ = Reclass_bare1
            rock_remap = f"1 {bare_rock_weight}"
            Reclass_bare1 = arcpy.sa.Reclassify(in_raster=x, reclass_field="VALUE", remap=rock_remap, missing_values="DATA")
            Reclass_bare1.save(Reclassify_3_)
            delete_list.append(Reclassify_3_)
        if "degraded" in x:
            # Process: Reclassify (Reclassify) (sa)
            Reclass_degraded2 = x
            Reclassify = Reclass_degraded2
            degraded_remap = f"1 {degraded_weight}"
            Reclass_degraded2 = arcpy.sa.Reclassify(in_raster=x, reclass_field="VALUE", remap=degraded_remap, missing_values="DATA")
            Reclass_degraded2.save(Reclassify)
            delete_list.append(Reclassify)
        if "crop" in x:
           # Process: Reclassify (2) (Reclassify) (sa)
           Reclass_crop5 = x
           Reclassify_2_ = Reclass_crop5
           crop_remap = f"1 {crop_weight}"
           Reclass_crop5 = arcpy.sa.Reclassify(in_raster=x, reclass_field="VALUE", remap=crop_remap, missing_values="DATA")
           Reclass_crop5.save(Reclassify_2_)
           delete_list.append(Reclassify_2_)
    # Process: Polygon to Raster (Polygon to Raster) (conversion)
    Mara_Agriculture_Clip_PolygonToRaster = maskpath+"_Raster"
    with arcpy.EnvManager(cellSize=cellsize):
        arcpy.conversion.PolygonToRaster(in_features=maskpath, value_field="OBJECTID",
                                         out_rasterdataset=Mara_Agriculture_Clip_PolygonToRaster, cell_assignment="CELL_CENTER",
                                         priority_field="NONE", cellsize=cellsize, build_rat="BUILD")
    arcpy.management.Delete(maskpath)
    # Process: Reclassify (4) (Reclassify) (sa)
    Reclass_Mara1 = Mara_Agriculture_Clip_PolygonToRaster
    Reclassify_4_ = Reclass_Mara1
    min = arcpy.sa.Raster(Mara_Agriculture_Clip_PolygonToRaster).minimum
    max = arcpy.sa.Raster(Mara_Agriculture_Clip_PolygonToRaster).maximum
    ag_remap = f"{min} {max} {agriculture_weight}"
    Reclass_Mara1 = arcpy.sa.Reclassify(in_raster=Mara_Agriculture_Clip_PolygonToRaster, reclass_field="VALUE",
                                        remap=ag_remap, missing_values="DATA")
    Reclass_Mara1.save(Reclassify_4_)
    delete_list.append(Reclassify_4_)
    # Mosaic to new Raster
    print("     Mosaic starting now")
    ag_combined_name1 = os.path.join(scratch, "LULC_Combo")
    ag_combined_name = "LULC_Combo"
    Mara_Ag_crop5_combined = arcpy.management.MosaicToNewRaster(input_rasters=[Reclass_Mara1, Reclass_crop5, Reclass_bare1, Reclass_degraded2],
                           output_location=scratch, raster_dataset_name_with_extension=ag_combined_name,
                           coordinate_system_for_the_raster= projection, pixel_type="8_BIT_UNSIGNED", cellsize=cellsize,
                           number_of_bands=1, mosaic_method="MAXIMUM", mosaic_colormap_mode="FIRST")[0]
    delete_list.append(ag_combined_name1)
    Mara_Ag_crop5_combined = arcpy.Raster(Mara_Ag_crop5_combined)
    with arcpy.EnvManager(mask=clip):
        arcpy.env.overwriteOutput = True
        # If scale is changed below has to be changed
        out_raster = arcpy.sa.Reclassify(Mara_Ag_crop5_combined, "Value", "1 1;2 2; 3 3; 4 4; 5 5;NODATA 0", "DATA")
        finalName = "LULC_FINAL"
        finaldata = os.path.join(output, finalName)
        LULC.name = finaldata
        out_raster.save(finaldata)
    delete_extras(delete_list)
    print(f"LCLU layer complete!\n")
    return finalName


def RoadsModel(primary, secondary, tertiary):
    print(f"Starting Roads Layer processing: ")
    primaryRoad = arcpy.sa.RemapRange(primary)
    secondaryRoad = arcpy.sa.RemapRange(secondary)
    tertiaryRoad = arcpy.sa.RemapRange(tertiary)
    remap_list = primaryRoad, secondaryRoad, tertiaryRoad

    # Check out necessary licenses
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("3D")

    delete_list = []
    extent = os.path.join(workspace, clip)
    arcpy.env.overwriteOutput = True

    # Process: Project (Project) (management)
    with arcpy.EnvManager(scratchWorkspace=scratch, workspace=workspace, extent=extent):
        print("     Projecting and Clipping layers")
        # Process: Project (Project) (management)
        out_data = polylines + "UTM"
        out_data = os.path.join(scratch, out_data)
        arcpy.management.Project(polylines, out_data, projection)
        delete_list.append(out_data)
        # Process: Clip (Clip) (analysis)
        out_data1 = polylines + "UTM_Clip"
        out_data1 = os.path.join(scratch, out_data1)
        arcpy.analysis.Clip(out_data, clip, out_data1)
        delete_list.append(out_data1)
        #extent = arcpy.Describe(clip).extent

        print("     Creating layers by road type")
        settlement_types = [['Road_Highway_Tarmac', 'Road_Primary_Dirt/Sand', 'Road_Primary_Murram', 'Road_Primary_Tarmac'],
                            ['Road_Secondary_Dirt/Sand', 'Road_Secondary_Murram', 'Road_Secondary_Tarmac', 'Road_Secondary_Unknown'],
                            ['Road_Tertiary_Dirt/Sand', 'Road_Tertiary_Murram', 'Road_Track_Dirt/Sand']]
        x = 1
        road_filename = []
        for settlement_type in settlement_types:
            filename = f"Road_{x}"
            filepath = os.path.join(scratch,filename)
            road_filename.append(filename)
            where = "', '".join(str(e) for e in settlement_type)
            x = x+1
            where_clause = f"type IN ('{where}')"
            arcpy.analysis.Select(in_features=out_data, out_feature_class=filepath,
                                  where_clause=f'{where_clause}')
            delete_list.append(filepath)

    with arcpy.EnvManager(workspace=scratch):
    #Run Distance Accumulation Tool on each road layer...
        print("     Running the Distance Accumulation tool on each road layer")
        # Process: Distance Accumulation (Distance Accumulation)
        distance_list = []
        for layer in road_filename:
            Distance_Accumulation = f"Distance_Acc_{layer}"
            with arcpy.EnvManager(cellSize=cellsize, extent=extent):
                Primary_Distance = arcpy.sa.DistanceAccumulation(in_source_data=layer, in_barrier_data="",
                                                             in_surface_raster="", in_cost_raster="", in_vertical_raster="",
                                                             vertical_factor="BINARY 1 -30 30", in_horizontal_raster="",
                                                             horizontal_factor="BINARY 1 45", out_back_direction_raster="",
                                                             out_source_direction_raster="",
                                                             out_source_location_raster="", source_initial_accumulation="",
                                                             source_maximum_accumulation="", source_cost_multiplier="",
                                                             source_direction="", distance_method="PLANAR")
                arcpy.management.Clip(in_raster=Primary_Distance, out_raster=Distance_Accumulation,
                                      in_template_dataset=os.path.join(workspace, clip),
                                      nodata_value="", clipping_geometry="ClippingGeometry",
                                      maintain_clipping_extent="MAINTAIN_EXTENT")
                #Primary_Distance.save(os.path.join(scratch, Distance_Accumulation))
                distance_list.append(Distance_Accumulation)
                delete_list.append(Distance_Accumulation)
    #Reclassify each road layer...

        print("     Reclassifying each road layer on a 0-5 scale")

    #Process: Reclassify (Reclassify) (sa)
        reclass_list = []
        for file, remap in zip(distance_list, remap_list):
            Reclassify = f"{file}_Reclass"
            PrimaryRC = arcpy.sa.Reclassify(in_raster=file, reclass_field="VALUE",
                                        remap=remap,
                                        missing_values="DATA")
            PrimaryRC.save(Reclassify)
            reclass_list.append(Reclassify)
            delete_list.append(Reclassify)

        print("     Mosaicking distance rasters")
        roads_combo = "roads_combo"
        RoadsCombined = arcpy.management.MosaicToNewRaster(input_rasters=reclass_list,
                           output_location=scratch, raster_dataset_name_with_extension=roads_combo,
                           coordinate_system_for_the_raster= projection, pixel_type="8_BIT_UNSIGNED", cellsize=cellsize,
                           number_of_bands=1, mosaic_method="MAXIMUM", mosaic_colormap_mode="FIRST")[0]
        delete_list.append(roads_combo)
        Road_Layer_Combo = arcpy.Raster(RoadsCombined)
    with arcpy.EnvManager(scratchWorkspace=scratch, workspace=workspace, mask=clip):
        arcpy.env.overwriteOutput = True
    # If scale is changed below has to be changed
        print("     Reclassifying mosaicked raster")
        out_raster = arcpy.sa.Reclassify(Road_Layer_Combo, "Value", "0 0;1 1;2 2; 3 3; 4 4; 5 5", "DATA")
        finalName = "Roads_FINAL"
        finaldata = os.path.join(output, finalName)
        RoadsModel.name = finaldata
        out_raster.save(finaldata)
    delete_extras(delete_list)
    print(f"Roads layer complete!\n")
    return finalName


def Albedo_Model_Mep():  # Albedo_Model_Mep
    print(f"Starting Albedo processing: ")
    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("spatial")
    delete_list = []

    print("     Projecting and Clipping layers")
    # Process: Clip Raster (Clip Raster) (management)
    albedo_clip = os.path.join(scratch, "Albedo_Clip")
    arcpy.management.Clip(in_raster=albedo, out_raster=albedo_clip, in_template_dataset=clip, nodata_value="1.79e+308",
                          clipping_geometry="NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
    albedo_clip = arcpy.Raster(albedo_clip)
    delete_list.append(albedo_clip)

    # Process: Project Raster (Project Raster) (management)
    albedo_project = os.path.join(scratch, "Albedo_Project")
    arcpy.management.ProjectRaster(in_raster=albedo_clip, out_raster=albedo_project, out_coor_system=projection,
                                           vertical="NO_VERTICAL")
    albedo_Project = arcpy.Raster(albedo_project)
    delete_list.append(albedo_Project)

    print("     Reclassifying albedo to 1-5 scale")
    # Process: Reclassify (Reclassify) (3d)
    albedo_reclass = os.path.join(scratch, "Albedo_Reclass")
    arcpy.ddd.Reclassify(in_raster=albedo_Project, reclass_field="VALUE",
                         remap="-17963.660573 -4707.755673 1;-4707.755673 235.124120 2;235.124120 6750.738393 3;"
                               "6750.738393 15063.763499 4;15063.763499 39328.809756 5", out_raster=albedo_reclass,
                         missing_values="NODATA")
    albedo_Reclass = arcpy.Raster(albedo_reclass)
    delete_list.append(albedo_Reclass)

    print(f"     Resampling albedo to cell size: {cellsize}")
    # Process: Resample (Resample) (management)
    finalName = "Albedo_FINAL"
    albedo_resample = os.path.join(output, "Albedo_FINAL")
    arcpy.management.Resample(in_raster=albedo_Reclass, out_raster=albedo_resample, cell_size=cellsize, resampling_type="NEAREST")
    albedo_resample = arcpy.Raster(albedo_resample)
    Albedo_Model_Mep.name = albedo_resample

    delete_extras(delete_list)
    print(f"Roads layer complete!\n")
    return finalName


def Pop_Density():  # Model 3
    print(f"Starting population density processing: ")
    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = workspace
    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("3D")

    pop_data = arcpy.Raster("Population_Data")
    delete_list = []
    # Process: Extract by Mask (Extract by Mask) (sa)
    pop_den_mask = os.path.join(scratch, "pop_den_mask")
    Extract_by_Mask = pop_den_mask
    pop_den_mask = arcpy.sa.ExtractByMask(in_raster=pop_data, in_mask_data=clip)
    pop_den_mask.save(Extract_by_Mask)
    delete_list.append(Extract_by_Mask)

    print("     Projecting and Clipping layers")
    # Process: Project Raster (Project Raster) (management)
    popden_project = os.path.join(scratch, "popden_project")
    arcpy.management.ProjectRaster(in_raster=pop_den_mask, out_raster=popden_project, out_coor_system=projection,
                                           vertical="NO_VERTICAL")
    popden_project = arcpy.Raster(popden_project)
    delete_list.append(popden_project)
    print(f"     Resampling layer to cell size: {cellsize}")
    # Process: Resample (Resample) (management)
    popden_resample = os.path.join(scratch, "popden_resample")
    arcpy.management.Resample(in_raster=popden_project, out_raster=popden_resample, cell_size=cellsize, resampling_type="NEAREST")
    popden_resample = arcpy.Raster(popden_resample)
    delete_list.append(popden_resample)

    print("     Reclassifying raster to 1-5 scale")
    # Process: Reclassify (Reclassify) (sa)
    finalName = "PopDen_FINAL"
    mara_reclass = os.path.join(output, finalName)
    Reclassify = mara_reclass
    Marareclass = arcpy.sa.Reclassify(in_raster=popden_resample, reclass_field="VALUE",
                                      remap="2.111691 12.976758 1;12.976758 13.978778 2;13.978778 24.843844 3;24.843844 142.655606 4;142.655606 1420.108521 5",
                                      missing_values="DATA")
    Marareclass.save(Reclassify)
    Pop_Density.name = Reclassify
    delete_extras(delete_list)
    return finalName


def Settlement(boma, town, safari, hotel):  # Model
    print(f"Starting settlement processing: ")

    boma_remap = arcpy.sa.RemapRange(boma)
    town_remap = arcpy.sa.RemapRange(town)
    safari_remap = arcpy.sa.RemapRange(safari)
    hotel_remap = arcpy.sa.RemapRange(hotel)
    remap_list = boma_remap, hotel_remap, town_remap, safari_remap

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = True

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")
    arcpy.CheckOutExtension("3D")

    arcpy.env.workspace = workspace

    landDx_points = points
    Eucl_distance_method = "PLANAR"
    delete_list = []

    print("     Projecting and Clipping layers")
        # Process: Clip (Clip) (analysis)
    landDx_points_Clip3 = os.path.join(scratch, "landDx_points_Clip")
    arcpy.analysis.Clip(in_features=landDx_points, clip_features=clip, out_feature_class=landDx_points_Clip3, cluster_tolerance="")
    delete_list.append(landDx_points_Clip3)

        # Process: Project (Project) (management)
    landDx_points_Project = os.path.join(scratch, "landDx_points_Project")
    arcpy.management.Project(in_dataset=landDx_points_Clip3, out_dataset=landDx_points_Project,
                             out_coor_system=projection, transform_method=[],
                             preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
    delete_list.append(landDx_points_Project)

    print(f'     Fixing data input: spaces and slashes are replaced with "_"')
    # Fixing data input so there are no spaces or /
    with arcpy.da.UpdateCursor(landDx_points_Project, "type") as cursor:
        for row in cursor:
            if " " in row[0]:
                name_fix = row[0].replace(" ", "_")
                row[0] = name_fix
                cursor.updateRow(row)
            if "/" in row[0]:
                name_fix = row[0].replace("/", "_")
                row[0] = name_fix
                cursor.updateRow(row)


    print("     Creating layers by settlement type")
    settlement_types = ['Settlement_Boma', 'Village', 'Settlement_Hotel_Lodge', 'Town', 'Settlement_Safari_Camp']
    settlement_types = ['Settlement_Boma', 'Settlement_Hotel_Lodge', 'Town', 'Settlement_Safari_Camp']
    for settlement_type in settlement_types:
        clausename = settlement_type
        typefilename = settlement_type.replace("/", "_")
        filename = os.path.join(scratch, typefilename)
        where_clause = f"type IN ('{clausename}')"
        arcpy.analysis.Select(in_features=landDx_points_Project, out_feature_class=filename,
                            where_clause=f'{where_clause}')
        delete_list.append(filename)

    print("     Running the distance accumulation tool on each settlement layer")
    reclass_list = []
    for filename, remap in zip(settlement_types, remap_list):
        EucDist = f"EucDist_{filename}"
        with arcpy.EnvManager(workspace=scratch, extent=os.path.join(workspace, clip)):
            EucDist_process = arcpy.sa.EucDistance(in_source_data=filename, cell_size=cellsize, distance_method=Eucl_distance_method)
            arcpy.management.Clip(in_raster=EucDist_process, out_raster=EucDist, in_template_dataset=os.path.join(workspace, clip),
                            nodata_value="", clipping_geometry="ClippingGeometry", maintain_clipping_extent="MAINTAIN_EXTENT")
            Reclassify_name = f"Reclass_{filename}"
            Reclassify = arcpy.sa.Reclassify(in_raster=EucDist, reclass_field="VALUE",
                                        remap=remap, missing_values="NODATA")
            Reclassify.save(Reclassify_name)
            delete_list.append(Reclassify_name)
            reclass_list.append(Reclassify_name)
            delete_list.append(EucDist)

    print("     Mosaicking distance rasters")
    # Process: Weighted Overlay (Weighted Overlay) (sa)
    settlementcombo = "settlementcombo"
    with arcpy.EnvManager(workspace=scratch, extent=os.path.join(workspace, clip)):
        SettlementsCombined = arcpy.management.MosaicToNewRaster(input_rasters=reclass_list, output_location=scratch,
                                                                 raster_dataset_name_with_extension=settlementcombo,
                                                                 coordinate_system_for_the_raster=projection,
                                                                 pixel_type="8_BIT_UNSIGNED", cellsize=cellsize,
                                                                 number_of_bands=1, mosaic_method="MAXIMUM",
                                                                 mosaic_colormap_mode="FIRST")[0]
        SettlementsCombined = arcpy.Raster(SettlementsCombined)
        delete_list.append(settlementcombo)
    # If overal scale is changed below has to be changed
        out_raster = arcpy.sa.Reclassify(SettlementsCombined, "Value", "0 0;1 1;2 2; 3 3; 4 4; 5 5", "DATA")
    # Saving Final Raster
        finaldataname = "Settlements_FINAL"
        finaldata = os.path.join(output, finaldataname)
        Settlement.name = finaldata
        out_raster.save(finaldata)
        delete_extras(delete_list)
        print(f"Settlement Layer complete\n")
        return finaldataname


def delete_extras(delete_list):
    with arcpy.EnvManager(workspace=scratch):
        if delete_extras_wanted:
            arcpy.management.Delete(delete_list)


def Weighted_Sum_Combination(fence_weight, LULC_weight, LandOwnership_weight, roads_weight, albedo_weight,
                 popDen_weight, settlement_weight):

    fence = "FenceWeight_FINAL"
    lulc = "LULC_FINAL"
    ownership = "LandOwnership_FINAL"
    roads = "Roads_FINAL"
    albedo = "Albedo_FINAL"
    settlement = "Settlements_FINAL"
    pop = "PopDen_FINAL"
    print("Summing all layers now:")
    arcpy.env.workspace = output
    arcpy.env.overwriteOutput = True
    out_raster = arcpy.sa.WeightedSum(arcpy.sa.WSTable([[fence, "Value", fence_weight],
                                                        [lulc, "Value", LULC_weight],
                                                        [ownership, "Value", LandOwnership_weight],
                                                        [roads, "Value", roads_weight],
                                                        [albedo, "Value", albedo_weight], [pop, "Value", popDen_weight],
                                                        [settlement, "Value", settlement_weight]]))
    HFI_FINAL = os.path.join(output, "HFI_FINAL")
    out_raster.save(HFI_FINAL)
    print(f"The final HFI can be found here: {HFI_FINAL}\n")


########################################################################################################################
''' Calling All Functions '''
def main():
    # FenceModel(wire_fence_weight=wire_fence_type, electric_fence_weight=electric_fence_type,
    #            other_fence_weight=other_fence_type)
    # LULC(bare_rock_weight=bare_rock_LCLU, crop_weight=crop_LCLU, degraded_weight=degraded_land_LCLU,
    #      agriculture_weight=agriculture_LCLU)
    # Land_Ownership(mr_weight=protected_land_type, cca_weight=community_conservancy_land_type,
    #                up_weight=unprotected_land_type)
    RoadsModel(primary=primary_road_type, secondary=secondary_road_type, tertiary=tertiary_road_type)
    # Albedo_Model_Mep()
    # Pop_Density()
    # Settlement(boma=boma_settlement_type, town=town_settlement_type, safari=safari_settlement_type,
    #            hotel=hotel_settlement_type)
    Weighted_Sum_Combination(fence_weight=Fence_Weight, LULC_weight=LCLU_Weight,
                             LandOwnership_weight=Land_Ownership_Weight,roads_weight=Roads_Weight,
                             albedo_weight=Albedo_Weight, popDen_weight=Population_Density_Weight,
                             settlement_weight=Settlement_Weight)
    print("Script Complete")


main()
