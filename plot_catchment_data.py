# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 13:30:06 2015

@author: sb708
"""

import cement
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import basename

from cement.core import foundation, controller
from cement.core.controller import expose
from cement.utils import shell

class GISAppController(controller.CementBaseController):
    class Meta:
        label = 'base'
        description = 'Python CLI application to automate QS data plots'
        arguments = [
            ( ['-d', '--data'], dict(action='store', dest='data',
                      help='path to data file') )]
        

class GISApp(foundation.CementApp):
    
    class Meta:
        label = 'QS_Plotter'
        base_controller = GISAppController


class Plotter():
    
    def __init__(self, datafile):
        self.precipitation = {}
        self.areas = {}
        self.areas_col = []
        self.relief = {}
        self.QS_t_yr = {}
        self.volume_m3_yr = {}
        self.erosion_mm_yr = {}
        self.erosion_col = []
        self.faults = {}
        self.distances = {}
        self.selected_fault = ''
        self.datalabel_x = '' 
        self.datalabel_y = ''
        self.dataunit_x = ''
        self.dataunit_y = ''
        self.plot_options = {
            'logarithmic': 0,
            'trend': 1,
            'line': 'None'
        }
        
        self.datanames = [
            'area',
            'precipitation',
            'relief',
            'volume',
            'sediment flux',
            'erosion m',
            'erosion mm'
        ]
        
        self.datalabels = [
            'Catchment areas',
            'Average precipitation',
            'Relief',
            'Sediment volume',
            'Sediment flux',
            'Catchment erosion',
            'Catchment erosion'
        ]
        
        self.dataunits = [
            'km$^2$',
            'mm/yr',
            'km',
            'm$^3$/kg',
            'T/yr',
            'm',
            'mm'
        ]
        
        with open(datafile, 'rb') as csvfile:
            fault_data = csv.reader(csvfile, delimiter=',')
            for row in fault_data:
                try:
                    if int(row[0]):
                        c_id = int(row[0])
                        self.precipitation.update({c_id: float(row[1])})
                        self.areas.update({c_id: float(row[5])})
                        self.areas_col.append(float(row[5]))
                        self.QS_t_yr.update({c_id: float(row[9])})
                        self.volume_m3_yr.update({c_id: float(row[11])})
                        self.erosion_mm_yr.update({c_id: float(row[13])})
                        self.erosion_col.append(float(row[13]))
                        self.distances.update({c_id: float(row[15])})                      
            
                        if int(row[14]) in self.faults:
                            self.faults[int(row[14])].append(c_id)
                        else:
                            self.faults.update({int(row[14]): [c_id]})
                
                except ValueError:
                    error = ValueError
        
        self.choose_plot()
        
    def choose_plot(self):
        
        x_axis_prompt = shell.Prompt("Choose X axis", options = self.datanames, numbered = True)
        x_choice = x_axis_prompt.input
        print(x_choice+' chosen for x-axis')
        y_axis_prompt = shell.Prompt("Choose Y axis", options = self.datanames, numbered = True)
        y_choice = y_axis_prompt.input
        print(y_choice+' chosen for y-axis')
        
        extent_choices = [
            'All catchments',
            'Fault specific',
            'Specific catchments'
        ]
        

        
        p1 = shell.Prompt("Data extent", options = extent_choices, numbered = True)
        
        if p1.input is 'Fault specific':
            fault_list = map(lambda x: 'fault '+str(x), self.faults.keys())
            p2 = shell.Prompt("Choose fault", options = fault_list, numbered = True)  
            
            self.selected_fault = int(p2.input.replace('fault ', ''))
            
            f_x_data, x_label_id, x_unit_id = self.get_fault_data(x_choice)
            f_y_data, y_label_id, y_unit_id = self.get_fault_data(y_choice)
            
            self.plot_data([f_x_data,x_label_id, x_unit_id], [f_y_data, y_label_id, y_unit_id])
            
        elif p1.input is 'Specific catchments':
            # This is when we plot everything
            print('yo')
            

    def select_dataset(self, dataname):
        print('Choosing dataset for '+ dataname)
            
        if dataname is self.datanames[1]:
            dataset = self.precipitation
            unit = 1
            label = 1
        elif dataname is self.datanames[0]:
            dataset = self.areas
            unit = 0
            label = 0
        elif dataname is self.datanames[2]:
            dataset = self.relief
            unit = 2
            label = 2
        elif dataname is self.datanames[3]:
            dataset = self.volume_m3_yr
            unit = 3
            label = 3
        elif dataname is self.datanames[4]:
            dataset = self.QS_t_yr
            self.plot_options['logarithmic'] = 1
            unit = 4
            label = 4
        elif dataname is self.datanames[6]:
            dataset = self.erosion_mm_yr
            label = 5
            unit = 5
        elif dataname is self.datanames[5]:
            dataset = self.erosion_m_yr
            unit = 6
            label = 6
        else:
            dataset = ''
            print(dataname+' - not recognised!')
            exit
            
        return dataset, unit, label
        

    def get_fault_data(self, dataname):
        
        dataset, label, unit = self.select_dataset(dataname)
        
        catchments = self.faults[self.selected_fault]
        
        output_data = []

        for c in catchments:
            output_data.append(dataset[c])

        return output_data, label, unit
    
    
    def plot_data(self, x, y):
        x_data = x[0]
        x_label = self.datalabels[x[1]]
        x_unit = self.dataunits[x[2]]
        
        y_data = y[0]
        y_label = self.datalabels[y[1]]
        y_unit = self.dataunits[y[2]]
        
        plt.plot(x_data, y_data, marker='o', linestyle=self.plot_options['line'])
        plt.xlabel(x_label+' ('+x_unit+')')
        plt.ylabel(y_label+' ('+y_unit+')')
        
        if self.plot_options['logarithmic']:
            plt.yscale('log')
            plt.xscale('log')
        

        if self.plot_options['trend']:
            
            if self.plot_options['logarithmic']:
                logx = np.log(x_data)
                logy = np.log(y_data)
                max_x = np.amax(logx)
                min_x = np.amin(logx)
                x_d = np.logspace(min_x, max_x)
                coeffs = np.polyfit(logx,logy,deg=1, full=True)
                print('Slope: '+str(coeffs[0][0]))
                print('Intercepts: '+str(coeffs[0][1]))
                poly = np.poly1d(coeffs[0])
                yfit = lambda x: np.exp(poly(np.log(x)))
                plt.loglog(x_d,yfit(x_d))
            else:
                
                m, b = np.polyfit(x_data, y_data, 1)
                x_min = np.amin(x_data)
                x_max = np.amax(x_data)
                x_d = np.linspace(np.floor(x_min), np.ceil(x_max))
                plt.plot(x_d, m*x_d + b, '-')
        
        plt.show()
        
app = GISApp()

app.setup()
    
app.run()

if app.pargs.data:
    if os.path.exists(app.pargs.data):
        plotter = Plotter(app.pargs.data)
