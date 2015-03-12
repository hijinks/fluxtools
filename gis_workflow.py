# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 15:14:57 2015

@author: Sam Brooke
"""
import cement
import yaml
import os
from os.path import basename

import arcpy
import datetime
import shutil
import math
import csv
import glob
from arcpy import env
from arcpy.sa import *

from cement.core import foundation, controller
from cement.core.controller import expose
from cement.utils import shell

class GISAppController(controller.CementBaseController):
    class Meta:
        label = 'base'
        description = 'Python CLI application to automate some GIS processing'
        arguments = [
            ( ['-c', '--config'], dict(action='store', dest='config',
                      help='path to config file') ),
            ( ['-b', '--batch'], dict(action='store', dest='batch',
                      help='path to batch directory') ),
            ]

    @expose(hide=True, aliases=['run'])
    def default(self):
        print("Running in default mode")

    @expose(help='Use existing batch and skip to watershed processing')
    def process_watersheds(self):
        print("Skipping to watersheds")
        self.app.skip_to_watersheds = 1

    @expose(help='Use existing batch and skip to watershed processing')
    def calculate_bqart(self):
        print("Skipping to discharge calculations")
        self.app.skip_to_watersheds = 1
        self.app.skip_to_discharge = 1        
        
class GISApp(foundation.CementApp):
    
    skip_to_watersheds = 0
    skip_to_discharge = 0
    
    class Meta:
        label = 'GIS_Automator'
        base_controller = GISAppController

   
app = GISApp()

class GISbatch:
    'Common base class for GIS batch processing'
   
    def __init__(self, config):
        self.project_root = config['root']
        self.project_name = config['project_name']
        self.projection_code = config['projection_code']
        self.pour_points_path = config['pour_points_path']
        self.scratch_path = config['scratch']
        self.output_path = config['output']
        self.original_dem = config['original_dem']
        
        # Workflow variables
        self.flow_dir = config['flow_dir']
        self.flow_acc = config['flow_acc']
        self.str_net = config['str_net']
        self.set_null = config['set_null']
        self.str_ord = config['str_ord']
        self.pour_points = config['pour_points']
        
        # Climate variables
        self.climates = config['climates']

        
        # Set the environment variables
        arcpy.env.scratchWorkspace = self.scratch_path
        self.batch_path = self.set_workspace()
        arcpy.env.workspace = self.batch_path
        
        sr = arcpy.SpatialReference(self.projection_code)
        arcpy.env.outputCoordinateSystem = sr

        # Load in Spatial Analyst Toolbox
        arcpy.CheckOutExtension("Spatial")
    
    def select_batch_directory(self, root_dir):
        os.chdir(root_dir)
        times = {}
        days = {}
        
        for dir_name in glob.glob("*"):
            t_frags = dir_name.split('_')
            t_int = map(int, t_frags)
            dt = datetime.datetime(t_int[0], t_int[1], t_int[2], t_int[3], t_int[4], t_int[5])
            times.update({dir_name: dt})
            ds = '_'.join([t_frags[0], t_frags[1], t_frags[2]])
            days.update({ds:{}})
                    
        for k in times.keys():
            t_frags = k.split('_')
            ds = '_'.join([t_frags[0], t_frags[1], t_frags[2]])
            t_string = times[k].strftime("%H:%M:%S")
            days[ds].update({t_string:k})
        
        return days
    
    def get_time_string(self):
        t = datetime.datetime.now()
        tstuff = [t.year, t.month, t.day, t.hour, t.minute, t.second]
        # Convert integer values to string
        tstring_list = map(str, tstuff)
        d = '_'.join(tstring_list)
        
        return d
        
    def set_workspace(self):
        
        if app.pargs.batch is None:
            
            dirname = self.get_time_string()
        
            # Assuming it doesn't exist already
            output_batch_path = os.path.join(self.output_path, dirname)
  
            os.makedirs(output_batch_path)
            
            # Copy original DEM
            shutil.copy2(self.original_dem, output_batch_path)
            
        else:
            if os.path.isdir(app.pargs.batch):
                output_batch_path = app.pargs.batch
            else:
                print('Batch directory does not exist')
                exit
        
        return output_batch_path 
        

    def hydro_workflow(self):
        print('Starting Hydrology Workflow...')
      
        print('Fill')
        fill_path = self.fill()
        
        print('Flow direction')
        flow_path = self.flow_direction(fill_path)
        
        print('Flow accumulation')
        flow_acc_path = self.flow_accumulation(flow_path)
        
        print('Steam network')
        stream_net_path = self.stream_network(flow_acc_path)

        print('Nullify')
        null_path = self.nullify(stream_net_path)
        
        print('Stream order')
        s_ord_path = self.stream_order(null_path, flow_path)
        
        print('Vectorise streams')
        vector_streams = self.vectorise_streams(s_ord_path, flow_path)
        
        # Save file values to YAML file
        hydro_paths = {
            'fill_path' : fill_path,
            'flow_path' : flow_path,
            'flow_acc_path' : flow_acc_path,
            'stream_net_path' : stream_net_path,
            'null_path' : null_path,
            's_ord_path' : s_ord_path,
            'vector_streams' : vector_streams
       }
        with open(os.path.join(self.batch_path,'hydro_paths.yml'), 'w') as outfile:
            outfile.write(yaml.dump(hydro_paths, default_flow_style=True) )
        
        return hydro_paths
        

    def watershed_workflow(self, original_pour_points, hydro_paths):

        print('Starting Watershed workflow')
      
        print('Creating batch directory')
        self.watershed_batch_path, pp_path = self.setup_watershed_batch(original_pour_points)

        print('Assigning UIDs') # For processing
        self.assign_uids(pp_path)
        
        print('Snap to pour points')
        snap_pp_path = self.snap_pour_points(pp_path, hydro_paths['flow_acc_path'])

        print('Extract watersheds')
        ws_path = self.watersheds(hydro_paths['flow_path'], snap_pp_path)

        print('Assigning CIDs') # Catchment Ids
        self.assign_cids(ws_path)
        
        # print('Converting to polygons')
        # self.ws_to_poly(ws_path)

        
        watershed_paths = {
            'pour_points' : snap_pp_path,
            'watersheds' : ws_path
        }
        
        with open(os.path.join(self.watershed_batch_path,'watershed_paths.yml'), 'w') as outfile:
            outfile.write(yaml.dump(watershed_paths, default_flow_style=True) )
            
        return ws_path


    def bqart_workflow(self, watershed_raster, hydro_paths, watershed_path, temp_directory, precip_directory):
        
        print('Creating climate batch directory')
        climate_batch_path = self.climate_batch_directory(watershed_path)
        originals_batch_path = os.path.join(climate_batch_path, 'originals')
        
        print('Averaging precipitation rasters')
        precip_raster_path = self.average_rasters(precip_directory, originals_batch_path, 'precip_combined.tif', 0)
        print('Averaging temperature rasters')        
        temp_raster_path = self.average_rasters(temp_directory, originals_batch_path, 'temp_combined.tif', 1)
        
        print('Clipping climate rasters') 
        temp_clip = self.clip_raster(temp_raster_path, originals_batch_path, 'temp_clip.tif', watershed_raster)
        precip_clip = self.clip_raster(precip_raster_path, originals_batch_path, 'precip_clip.tif', watershed_raster)
        
        print('Climate zone statistics') 
        tz_dat_path = self.zone_statistics(climate_batch_path, watershed_raster, temp_clip, 'temp_data')
        pz_dat_path = self.zone_statistics(climate_batch_path, watershed_raster, precip_clip, 'precip_data')
        ez_dat_path = self.zone_statistics(climate_batch_path, watershed_raster, hydro_paths['fill_path'], 'elev_data')
        
        print('Calculating Qs using BQART') 
        qs_data = self.do_bqart(pz_dat_path, tz_dat_path, ez_dat_path)
        
        self.save_data_to_csv(qs_data, climate_batch_path)
        
    # ARC GIS PROCESSES
    # Hydro stuff

    def fill(self):
        fill_z_limit = ""

        out_fill = Fill(self.original_dem, fill_z_limit)
        out_fill_raster = self.project_name + '_fill.tif'
        out_fill_path = os.path.join(self.batch_path, out_fill_raster)
        out_fill.save(out_fill_path)
        
        return out_fill_path
        

    def flow_direction(self, fill_path):
        force_flow = self.flow_dir['force_flow']
        
        out_flow_dir = FlowDirection(fill_path, force_flow)
        out_flow_dir_raster = self.project_name + '_f_dir.tif'
        out_flow_dir_path = os.path.join(self.batch_path, out_flow_dir_raster)
        out_flow_dir.save(out_flow_dir_path)
        
        return out_flow_dir_path
        

    def flow_accumulation(self, flow_path):
        flow_weight_raster = self.flow_acc['flow_weight_raster']
        flow_data_type = self.flow_acc['flow_data_type']

        out_flow_acc = FlowAccumulation(flow_path, flow_weight_raster, flow_data_type)
        out_flow_acc_raster = self.project_name + '_f_acc.tif'
        out_flow_acc_path = os.path.join(self.batch_path, out_flow_acc_raster)
        out_flow_acc.save(out_flow_acc_path)
        
        return out_flow_acc_path
        

    def stream_network(self, flow_acc_path):
        con_where_clause = self.str_net['conditional']
        false_constant = self.str_net['false_constant']
        true_raster = flow_acc_path

        out_con = Con(flow_acc_path, true_raster, false_constant, con_where_clause)
        out_con_raster = self.project_name + '_net.tif'
        stream_net_path = os.path.join(self.batch_path, out_con_raster)
        out_con.save(stream_net_path)
        
        return stream_net_path
 
    def nullify(self, stream_net_path):
        false_raster = self.set_null['false_raster']
        null_where_clause = self.set_null['conditional']

        out_null = SetNull(stream_net_path, false_raster, null_where_clause)
        out_null_raster = self.project_name + '_net_null.tif'
        out_null_path = os.path.join(self.batch_path, out_null_raster)
        out_null.save(out_null_path)
        
        return out_null_path       

    def stream_order(self, null_path, flow_path):
        method = self.str_ord['method']
        
        out_s_ord_raster = self.project_name + '_s_order.tif'
        out_s_ord_path = os.path.join(self.batch_path, out_s_ord_raster)
        out_s_ord = StreamOrder(null_path, flow_path, method)
        out_s_ord.save(out_s_ord_path)
        
        return out_s_ord_path
        

    def vectorise_streams(self, s_ord_path, flow_path):
        out_sf_name = self.project_name + '_streams'
        out_sf_path = os.path.join(self.batch_path, out_sf_name)
        StreamToFeature(s_ord_path, flow_path, out_sf_path)
        
        return out_sf_path
        
    
    
    # Watershed stuff
    
    def setup_watershed_batch(self, original_pour_points):
        # Each watershed calculations need to be discrete from one another
        timestamp = self.get_time_string()
        
        if os.path.isdir(os.path.join(self.batch_path, 'watershed_calcs')) == 0:
            os.makedirs(os.path.join(self.batch_path, 'watershed_calcs'))
            
        watershed_batch_path = os.path.join(self.batch_path, 'watershed_calcs', timestamp)
        originals_batch_path = os.path.join(watershed_batch_path, 'originals')        
        os.makedirs(watershed_batch_path)
        os.makedirs(originals_batch_path)
        
        # Copy original Pour Points
        shutil.copy2(original_pour_points, originals_batch_path)
        
        working_pp_path = os.path.join(watershed_batch_path, basename(original_pour_points))
        arcpy.CopyFeatures_management(original_pour_points, working_pp_path)
        
        return watershed_batch_path, working_pp_path
            
            
    def assign_uids(self, pp_path):
        arcpy.AddField_management(pp_path, 'UID', "SHORT")
        rows = arcpy.da.UpdateCursor(pp_path, ['UID'])
        i = 1
        for row in rows:
            row[0] = i
            i = i+1
            rows.updateRow(row)
       
       # Unlock data
        del row 
        del rows
 
    def assign_cids(self, pp_path):
        arcpy.AddField_management(pp_path, 'c_id', "TEXT")
        rows = arcpy.da.UpdateCursor(pp_path, ['c_id'])
        i = 1
        for row in rows:
            row[0] = 'c_'+str(i)
            i = i+1
            rows.updateRow(row)
       
       # Unlock data
        del row 
        del rows
        
    def pour_points_to_raster(self, pour_points):
        pp_raster_name = self.project_name +'_pp_raster.tif'
        pp_raster_path = os.path.join(self.watershed_batch_path, pp_raster_name)
        arcpy.PointToRaster_conversion(pour_points, "UID", pp_raster_path, 'MOST_FREQUENT', '', '10')        
        
        return pp_raster_path
        
        
    def snap_pour_points(self, pour_points, flow_acc):
        snap_distance = self.pour_points['snap_distance']
        
        out_pp_name = self.project_name + '_snap_ppoints.tif'
        out_pp_path = os.path.join(self.watershed_batch_path, out_pp_name)        
        pp = SnapPourPoint(pour_points, flow_acc, snap_distance, "UID")
        pp.save(out_pp_path)
        
        return out_pp_path
    
            
    def watersheds(self, flow_path, pp_path):
        inPourPointField = "VALUE" # Now contains the c_id values
        
        out_ws_name = self.project_name + '_watersheds.tif'
        out_ws_path = os.path.join(self.watershed_batch_path, out_ws_name)
        outWatershed = Watershed(flow_path, pp_path, inPourPointField)
        outWatershed.save(out_ws_path)
        
        return out_ws_path
        
        
    def ws_to_poly(self, ws_path):
        
        out_poly_name = self.project_name + '_poly_ws'
        out_poly_path = os.path.join(self.watershed_batch_path, out_poly_name)
        arcpy.RasterToPolygon_conversion(ws_path, out_poly_path, "NO_SIMPLIFY", 'VALUE')      
        
        return out_poly_path
        
    
    # BQART stuff
    
        
    def climate_batch_directory(self, watershed_directory):
        timestamp = self.get_time_string()
        print('Creating batch files')
        climate_batch_path = os.path.join(watershed_directory, 'climate_calcs', str(timestamp))
        
        os.makedirs(climate_batch_path)
        climate_originals_batch_path = os.path.join(climate_batch_path, 'originals')  
        os.makedirs(climate_originals_batch_path)
        
        return climate_batch_path


    def average_rasters(self, search_directory, save_directory, name, monthly):
        os.chdir(search_directory)
        rasters = []
        for file in glob.glob("*.tif"):
            rasters.append(Raster(os.path.join(search_directory,file)))
         
        raster_sum = sum(rasters)
        
        if monthly: # Temp
            n_rasters = len(rasters)
            combined_raster = raster_sum / n_rasters
        else: # precip
            combined_raster = raster_sum

        combined_raster_path = os.path.join(save_directory, name)
        combined_raster.save(combined_raster_path)
        
        return combined_raster_path


    def clip_raster(self, input_raster, save_directory, name, extent):
        clip_raster_path = os.path.join(save_directory, name)
        arcpy.Clip_management(input_raster, '#', clip_raster_path, extent)    
        
        return clip_raster_path
    

    def zone_statistics(self, table_directory, watersheds, value_raster, data_name):
        table_path = os.path.join(table_directory, data_name)
        outdata = ZonalStatisticsAsTable(watersheds, "c_id", value_raster, table_path, "DATA")
        
        return outdata      
    
        
    def do_bqart(self, pz_data, tz_data, ez_data):
        
        t_cursor = arcpy.SearchCursor(tz_data)
        p_cursor = arcpy.SearchCursor(pz_data)
        e_cursor = arcpy.SearchCursor(ez_data)
        
        temps = {}
        precips = {}
        max_reliefs = {}
        min_reliefs = {}
        areas = {}
        
        for row in t_cursor:
            # Get mean temperature
            temps.update({row.getValue('c_id'): row.getValue('MEAN')})
        
        for row in p_cursor:
            # Get mean precipitation
            precips.update({row.getValue('c_id'): row.getValue('MEAN')})
            
        for row in e_cursor:
            # Get highest, lowest elevation & area of waters====heds
            max_reliefs.update({row.getValue('c_id'): row.getValue('MAX')})
            min_reliefs.update({row.getValue('c_id'): row.getValue('MIN')})
            areas.update({row.getValue('c_id'): row.getValue('AREA')})
            
        # BQART
        w = 0.0006
        B = 1
        
        qs_rows = []

        
        # Units!!

        # Qs (kg/s)
        # Qs (m^3/s)
        # A (km^2)
        # R (km)
        # T (C)
        
        # precips are in mm
        # temps are in C x 10
        # relief is in m
        # area is m^2
        
        
        for k in precips.keys():

            # Multiply area (m^2) with mean annual precip (m) 
            Q = math.pow(areas[k]*(precips[k]/float(1000)), 0.31)
            # Converting area from m^2 to km^2
            A = math.pow(areas[k]/float(1000000), 0.5)
            # Converting elevation from m to km
            R = max_reliefs[k] - min_reliefs[k]
            R_km = R/ float(1000)
            T = temps[k]/10 # Worldclim temps need to be divided by 10
            Qs = w*B*Q*A*R_km*T
            qs = [k, w, B, Q, A, R_km, T, Qs]
            qs_rows.append(qs)
            
        return qs_rows
    
    def save_data_to_csv(self, qs_data, path):
        data_name = 'qs_data.csv'
        row_headers = ['id', 'w', 'B', 'Q (kg/s)', 'A (km^2)', 'R (km)', 'T(C)', 'Qs (MT/y)']
        data_path = os.path.join(path, data_name)
        with open(data_path, 'wb') as qs_file:
            a = csv.writer(qs_file, delimiter=',')
            a.writerow(row_headers)
            for r in qs_data:
                a.writerow(r)
         
        print('Data saved to '+data_path)
        
try:
    app.setup()
    
    app.run()
    
    if app.pargs.config:
        try:
            f = open(app.pargs.config)
            yaml_config = yaml.load(f.read())
            f.close()
            gbatch = GISbatch(yaml_config)
            
            if app.skip_to_discharge == 0:
                if app.skip_to_watersheds == 0:
                    hydro_paths = gbatch.hydro_workflow()
                else:
                    try:
                        hydro_file_path = os.path.join(gbatch.batch_path, 'hydro_paths.yml')
                        
                        if os.path.exists(hydro_file_path):
                            f = open(hydro_file_path)
                            hydro_paths = yaml.load(f.read())
                            f.close()
                        else:
                            print('Cannot find '+ hydro_file_path)
                            hydro_paths_exists = 0
                            while hydro_paths_exists == 0:
                                p = shell.Prompt("Path to hydro_paths config file: ")
                                if os.path.exists(p.input):
                                    f = open(p.input)
                                    hydro_paths = yaml.load(f.read())
                                    f.close() 
                                    hydro_paths_exists = 1
                                else:
                                    print('File does not exist!')
                        
                    except (OSError, IOError) as e:
                        print(e)
                        exit

                if os.path.exists(gbatch.pour_points_path):
                    pour_point_path = gbatch.pour_points_path
                else:
                    pour_point_path = 0
                    while pour_point_path == 0:
                        p = shell.Prompt("Path to pour point shapefile: ")
                        if os.path.exists(p.input):
                            pour_point_path = p.input
                        else:
                            print('File does not exist!')

                watershed_raster = gbatch.watershed_workflow(pour_point_path, hydro_paths)
                watershed_directory = os.path.dirname(os.path.realpath(watershed_raster))
                
            else: # Skip to discharge calculations

                hydro_paths = 0
                watershed_raster = 0
                
                while hydro_paths == 0:
                    p = shell.Prompt("Path to hydro_paths config file: ")
                    if os.path.exists(p.input):
                        h_paths = p.input
                        hydro_paths = 1
                    else:
                        print('File does not exist!')
                f = open(h_paths)
                hydro_paths = yaml.load(f.read())
                f.close()   
                
                h_dir = os.path.dirname(os.path.realpath(h_paths))
                watershed_calcs = os.path.join(h_dir, 'watershed_calcs')
                days = gbatch.select_batch_directory(watershed_calcs)
                
                # Choose watershed batch
                d_strings = []
                d_frag_list = []
                
                # Ask which days
                for k in days:
                    d_frags = k.split('_')
                    d_int = map(int, d_frags)
                    dt = datetime.datetime(d_int[0], d_int[1], d_int[2])
                    d_strings.append(dt.strftime("%a %b %d %Y"))
                    d_frag_list.append(k)
                
                
                p1 = shell.Prompt("Pick batch day", options = d_strings, numbered = True)
                    
                # Ask which times
                d_index = d_strings.index(p1.input)
                d_choice = days[d_frag_list[d_index]]
                                    
                p2 = shell.Prompt("Pick batch time", options = d_choice.keys(), numbered = True)
                
                watershed_raster = os.path.join(watershed_calcs, 
                                                d_choice[p2.input], 
                            gbatch.project_name + '_watersheds.tif')
                
                watershed_directory = os.path.dirname(os.path.realpath(watershed_raster))
                    
            # Pick climate scenario
            climate_by_name = {}
            climate_names = []
            for c in gbatch.climates:
                climate_by_name.update({c['name']: c})
                climate_names.append(c['name'])
            
            p = shell.Prompt("Pick climate scenario", options = climate_names, numbered = True)

            climate_scenario = climate_by_name[p.input]
            temperature_directory = climate_scenario['temp_directory']
            precipitation_directory = climate_scenario['precip_directory']
            
            gbatch.bqart_workflow(watershed_raster, hydro_paths, 
                                  watershed_directory, temperature_directory, 
                                  precipitation_directory)
            
        except (OSError, IOError) as e:
            print(e)
            exit
    else:
        print('Please define path to config file -c CONFIG')
    
finally:
    arcpy.CheckInExtension("Spatial")
    app.close()
 


        
        