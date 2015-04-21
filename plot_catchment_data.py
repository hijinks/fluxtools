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
        self.erosion_m_yr = {}
        self.erosion_col = []
        self.faults = {}
        self.distances = {}
        self.selected_fault = ''
        self.datalabel_x = '' 
        self.datalabel_y = ''
        self.dataunit_x = ''
        self.dataunit_y = ''
        self.title = ''
        self.plot_options = {
            'logarithmic': { 'x': False, 'y': False },
            'trend': 1,
            'line': 'None'
        }
        
        self.datanames = [
            'area',
            'precipitation',
            'relief',
            'volume',
            'sediment flux',
            'erosion mm',
            'erosion m',
            'fault distance'
        ]
        
        self.datalabels = [
            'Catchment areas',
            'Average precipitation',
            'Relief',
            'Sediment volume',
            'Sediment flux',
            'Catchment erosion',
            'Catchment erosion',
            'Distance along fault'
        ]
        
        self.dataunits = [
            'km$^2$',
            'mm/yr',
            'km',
            'm$^3$/kg',
            'T/yr',
            'mm',
            'm',
            'm'
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
                        self.erosion_m_yr.update({c_id: float(row[12])})
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
            self.title = 'Fault '+str(self.selected_fault)
            
            f_x_data, x_label_id, x_unit_id = self.get_fault_data(x_choice, 'x')
            f_y_data, y_label_id, y_unit_id = self.get_fault_data(y_choice, 'y')
            
            s_x_data, s_y_data = self.sort_by_x(f_x_data, f_y_data)
            
            self.plot_data([s_x_data,x_label_id, x_unit_id], [s_y_data, y_label_id, y_unit_id])
            
        elif p1.input is 'Specific catchments':
            # This is when we plot everything
            print('yo')
        else:
            print('All catchments!')
            self.title = 'All Catchments'
            r_x_data, x_label_id, x_unit_id = self.select_dataset(x_choice, 'x')
            r_y_data, y_label_id, y_unit_id = self.select_dataset(y_choice, 'y')
            f_x_data = r_x_data.values()
            f_y_data = r_y_data.values()
            s_x_data, s_y_data = self.sort_by_x(f_x_data, f_y_data)
            self.plot_data([s_x_data,x_label_id, x_unit_id], [s_y_data, y_label_id, y_unit_id])

            
    def sort_by_x(self,xdata,ydata):
        unsorted = []
        
        def getKey(item):
            return item[0]
            
        for n in range(0, len(xdata)):
            unsorted.append([xdata[n], ydata[n]])
        
        sorted_by_x = sorted(unsorted, key=getKey)
        
        x_data_sorted = []
        y_data_sorted = []
        
        for m in range(0, len(sorted_by_x)):
            x_data_sorted.append(sorted_by_x[m][0])
            y_data_sorted.append(sorted_by_x[m][1])
            
        return x_data_sorted, y_data_sorted
        
        
    def select_dataset(self, dataname,x_y):
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
            
            if self.plot_options['logarithmic']['x'] is not 'force_off':
                self.plot_options['logarithmic']['x'] = True
                
            if self.plot_options['logarithmic']['y'] is not 'force_off':
                self.plot_options['logarithmic']['y'] = True
                
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
        elif dataname is self.datanames[7]:
            dataset = self.distances
            self.plot_options['logarithmic'][x_y] = 'force_off'
            self.plot_options['trend'] = 0
            self.plot_options['line'] = '-'
            unit = 7
            label = 7
        else:
            dataset = ''
            print(dataname+' - not recognised!')
            exit
            
        return dataset, unit, label
        

    def get_fault_data(self, dataname, x_y):
        
        dataset, label, unit = self.select_dataset(dataname, x_y)
        
        catchments = self.faults[self.selected_fault]
        
        output_data = []

        for c in catchments:
            output_data.append(dataset[c])

        return output_data, label, unit
    
    def get_r_squared(self, x_data, y_data, m, b):
        
        y_mean = sum(y_data)/len(y_data)        
        
        ss_totals = []
        
        for y in y_data:
            ss_totals.append(np.square(y-y_mean))
        
        ss_total = sum(ss_totals)
        
        f = []
        
        for x in x_data:
            f.append((m*x + b))

        ss_residuals = []
        
        
        y_f = np.array([y_data,f])
        
        
        for k in range(0,len(y_data)):
            ss_residuals.append(np.square(y_f[0][k]-y_f[1][k]))
        
        ss_res = sum(ss_residuals)
        
        r_squared = 1 - (ss_res/ss_total)
        
        return r_squared
            
            
    def plot_data(self, x, y):
        x_data = x[0]
        x_label = self.datalabels[x[1]]
        x_unit = self.dataunits[x[2]]
        
        y_data = y[0]
        y_label = self.datalabels[y[1]]
        y_unit = self.dataunits[y[2]]
        
        plt.plot(x_data, y_data, marker='o', linestyle=self.plot_options['line'])
        plt.title(self.title)
        plt.xlabel(x_label+' ('+x_unit+')')
        plt.ylabel(y_label+' ('+y_unit+')')
        
        log_settings = self.plot_options['logarithmic']
        
        loglog = False
        log_x = False
        log_y = False

        if log_settings['x'] is True:
            log_x = True
            
        if log_settings['y'] is True:
            log_y = True
            
        if log_settings['x'] is 'force_off':
            loglog = False
            log_x = False

        if log_settings['y'] is 'force_off':
            loglog = False
            log_y = False
        
        if log_x and log_y:
            loglog = True
            
        if loglog:
            plt.yscale('log')
            plt.xscale('log')
        
        if log_x:
            plt.xscale('log')
        
        if log_y:
            plt.yscale('log')
        
        
        if self.plot_options['trend']:
            
            if loglog:
                logx = np.log(x_data)
                logy = np.log(y_data)
                max_x = np.amax(logx)
                max_y = np.amax(logy)
                min_x = np.amin(logx)
                x_d = np.logspace(min_x, max_x)
                coeffs = np.polyfit(logx,logy,deg=1, full=True)
                print('Slope: '+str(coeffs[0][0]))
                print('Intercepts: '+str(coeffs[0][1]))
                poly = np.poly1d(coeffs[0])
                yfit = lambda x: np.exp(poly(np.log(x)))
                plt.loglog(x_d,yfit(x_d))
                r_squared = self.get_r_squared(logx,logy, coeffs[0][0], coeffs[0][1])
                slope = coeffs[0][0]
                intercept = coeffs[0][1]
                anno_x = max_x/3
                anno_x_t = anno_x*+(max_x*0.2)
                anno_y = max_y/3
                anno_y_t = anno_y+(max_y*0.2)
                print('r squared '+str(r_squared))
            elif log_x:
                print('Log x')
            elif log_y:
                print('Log y')
            else:
                m, b = np.polyfit(x_data, y_data, 1)
                r_squared = self.get_r_squared(x_data,y_data, m, b)
                print('r squared '+str(r_squared))
                x_min = np.amin(x_data)
                x_max = np.amax(x_data)
                y_max = np.amax(y_data)
                slope = m
                intercept = b
                anno_x = x_max/3
                anno_x_t = anno_x+(x_max*0.2)
                anno_y = y_max/3
                anno_y_t = anno_y+(y_max*0.2)
                x_d = np.linspace(np.floor(x_min), np.ceil(x_max))
                plt.plot(x_d, m*x_d + b, '-')
                
            if r_squared:
                plt.annotate('R$^2$ = '+str(round(r_squared, 3))+'\nSlope ='+str(slope)+
                '\nIntercept = '+str(intercept),  
                     xy=(anno_x, anno_y), xycoords='data',
                     xytext=(anno_x_t,anno_y_t), textcoords='data')
                 
        plt.show()
        
        
app = GISApp()

app.setup()
    
app.run()

if app.pargs.data:
    if os.path.exists(app.pargs.data):
        plotter = Plotter(app.pargs.data)
