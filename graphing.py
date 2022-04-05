# -*- coding: utf-8 -*-
"""

Graphing class design to make graph objects.
Input can be dataframe or arrays.

- set colours and themes for ions and voltages
- set style of no major gridlines etc.
- integration with bqplot for user interface


@author: E Shorer
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

dict_color ={"cl":'#1b9e77',"k":'#7570b3','x':'#d95f02','na':'#e7298a' }

clcolor='#1b9e77'
kcolor='#7570b3'
xcolor='#d95f02'
nacolor='#e7298a'
wcolor='k'



class graph(object):

    def __init__(self, _title = '', num_bins=200, start_time=45):
        """
        Basic settings of the graphing module
        Start time is the time to start plotting in seconds
        """
        self.title = _title
        self.bins = num_bins
        self.start_t = start_time*1000

    def graph_time_vm_allcomps(self, _vm_dict=dict,_time_arr=[]):
        """
        Displays the voltages for all the compartments versus time on a line graph
        _vm_arr is a dictionary of all the compartment voltages:
            {"comp name":[vm_array], ... }
        _time_arr is a list of all the time values
        """
        t_arr = _time_arr
        count = len(_time_arr)-self.start_t
        self.vm_dict = _vm_dict
        self.num_comps = len(_vm_dict)

        fig,ax = plt.subplots()
        ax.set_title("Membrane potential vs. Time per compartment")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Membrane Potential (mV)")

        x_arr =[]
        y_arr = []


        for i in range(self.bins):
            x_arr.append(t_arr[self.start_t+(round(count/self.bins)*i)])

        for k,v in self.vm_dict.items():
            temp_arr =[]
            comp_name =k
            y_arr =[]
            temp_arr = v
            for b in range(self.bins):
                y_arr.append(temp_arr[self.start_t+(round(count/self.bins)*b)])
            ax.plot(x_arr,y_arr,label=comp_name)


        ax.legend()

        return [fig,ax]




    def final_vals(self, dataframe):
        self.df = dataframe

    def set_x_axis(self, x_arr, x_title="x axis",):
        self.x_arr = x_arr
        self.x_title = x_title

    def set_y_axis(self, y_arr, y_title="y axis"):
        self.y_arr = y_arr
        self.y_title = y_title

    def plot_it(self,share_x ="True",share_y ="True"):

        plt.plot(self.x_arr, self.y_arr,)
        sns.despine()

