#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 19:39:59 2019

@author: doru
"""

from os.path import join, exists
from os import listdir, mkdir
from shutil import rmtree
import cv2
import numpy as np
import tkinter as tk
from PIL import Image, ImageTk
import pandas as pd

class GeneViewer(object):
    def __init__(self, feature_addrs, clustering_addr, gene_info_addr):
        self.feature_addrs = feature_addrs
        self.clustering_addr = clustering_addr
        self.gene_info_addr = gene_info_addr
        
        self.clustering = pd.read_csv(self.clustering_addr)
        self.gene_info  = pd.read_csv(self.gene_info_addr)
        
        self.clusters = self.clustering.Cluster.unique()
        np.ndarray.sort(self.clusters)
        
        self.data = {}
        for cluster in self.clusters:
            gene_names = [gene_name for gene_name in self.clustering.GeneNames[self.clustering.Cluster == cluster]]
            cluster_name = "Cluster_{number}".format(number = cluster)
            if not exists(join('group_descriptions', cluster_name)):
                fobj = open(join('group_descriptions', cluster_name), 'w')
                cluster_description = ""
                fobj.writelines(cluster_description)
            else:
                fobj = open(join('group_descriptions', cluster_name), 'r')
                cluster_description = fobj.read().strip()
            fobj.close()
            self.data[cluster_name] = [gene_names, cluster_description]
        
        self.current_cluster_index = 0 
        self.current_cluster = int(self.clusters[self.current_cluster_index])
        self.current_gene_index = 0 
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        
        self.root = tk.Tk()
        
        self.groups_frame      = tk.Frame(self.root)
        self.genes_frame       = tk.Frame(self.root)
        self.panel_frame       = tk.Frame(self.root)
        self.description_frame = tk.Frame(self.root)
        
        self.groups_frame.grid(     row = 0, column = 0)
        self.genes_frame.grid(      row = 0, column = 1)
        self.description_frame.grid(row = 0, column = 2, sticky = tk.N)
        self.panel_frame.grid(      row = 0, column = 3)
        
        self.buttons_frame   = tk.Frame(self.panel_frame)
        self.gene_name_label = tk.Label(self.panel_frame)
        self.canvas          = tk.Canvas(self.panel_frame, width = 500, height = 500)
        
        self.buttons_frame.grid(row = 0, column = 0)
        self.gene_name_label.grid(row = 1, column = 0)
        self.canvas.grid(row = 2, column = 0)
        
        self.save_button                  = tk.Button(self.buttons_frame, text = 'Save changes')
        self.changeGroupAssignment_button = tk.Button(self.buttons_frame, text = 'Change group assignment')
        self.merge_groups_button          = tk.Button(self.buttons_frame, text = 'Merge with')
        
        self.save_button.grid(                 row = 0, column = 0)
        self.changeGroupAssignment_button.grid(row = 0, column = 1)
        self.merge_groups_button.grid(         row = 0, column = 2)
        
        self.group_label = tk.Label(self.groups_frame)
        self.group_list_frame = tk.Frame(self.groups_frame)
        self.group_label.grid(     row = 0, column = 0)
        self.group_list_frame.grid(row = 1, column = 0)
        
        self.group_list = tk.Listbox(self.group_list_frame, height = 30, exportselection = 0)
        self.group_list.config(selectmode = tk.SINGLE)
        self.group_list_scroll = tk.Scrollbar(self.group_list_frame)
        
        self.group_list_scroll.pack(side = tk.RIGHT, fill = tk.Y)
        self.group_list.pack()
        self.group_list.config(yscrollcommand = self.group_list_scroll.set)
        self.group_list_scroll.config(command = self.group_list.yview)
        
        self.group_label['text'] = "Gene groups"
        
        self.gene_label = tk.Label(self.genes_frame)
        self.gene_list_frame = tk.Frame(self.genes_frame)
        self.gene_label.grid(     row = 0, column = 0)
        self.gene_list_frame.grid(row = 1, column = 0)
        
        self.gene_list = tk.Listbox(self.gene_list_frame, height = 30, exportselection = 0)
        self.gene_list.config(selectmode = tk.SINGLE)
        self.gene_list_scroll = tk.Scrollbar(self.gene_list_frame)
        
        self.gene_list_scroll.pack(side = tk.RIGHT, fill = tk.Y)
        self.gene_list.pack()
        self.gene_list.config(yscrollcommand = self.gene_list_scroll.set)
        self.gene_list_scroll.config(command = self.gene_list.yview)
        
        self.query_frame = tk.Frame(self.description_frame)
        self.group_description_label   = tk.Label(self.description_frame)
        self.group_description_content = tk.Text(self.description_frame, width = 50, height = 10, borderwidth=2, relief="solid")
        self.gene_description_label    = tk.Label(self.description_frame)
        self.gene_description_content  = tk.Text(self.description_frame, width = 50, height = 20, borderwidth=2, relief="solid")
        
        self.query_frame.grid(              row = 0, column = 0, sticky = tk.W)
        self.group_description_label.grid(  row = 1, column = 0, sticky = tk.W)
        self.group_description_content.grid(row = 2, column = 0, sticky = tk.W)
        self.gene_description_label.grid(   row = 3, column = 0, sticky = tk.W)
        self.gene_description_content.grid( row = 4, column = 0, sticky = tk.W)
        
        self.query_label = tk.Label(self.query_frame)
        self.query_enter = tk.Entry(self.query_frame, width = 15)
        
        self.query_label.grid(row = 0, column = 0)
        self.query_enter.grid(row = 0, column = 1)
        self.query_label['text'] = 'Enter gene name: '
        
        self.group_description_label['text'] = "Group description: "
        self.gene_description_label['text']  = "Gene summary: "
        
        self.load_image()
        
        self.update_genes()
        self.update_groups()
        
        self.gene_list.bind(                   "<ButtonRelease-1>", self.select_gene)
        self.group_list.bind(                  "<ButtonRelease-1>", self.select_group)
        self.save_button.bind(                 "<ButtonRelease-1>", self.save_changes)
        self.changeGroupAssignment_button.bind("<ButtonRelease-1>", self.changeGroupAssignment)
        self.merge_groups_button.bind(         "<ButtonRelease-1>", self.merge_groups)
        self.query_enter.bind(                 "<Return>",          self.query_gene_name)
        
        self.root.bind("<Left>",  self.go_to_previous)
        self.root.bind("<Right>", self.go_to_next)
        
        self.root.bind("<Up>",   self.previous_group)
        self.root.bind("<Down>", self.next_group)
        
        self.group_description_content.bind("<KeyRelease>", self.update_group_description)
        
        self.root.title("Gene expression with grouped genes")
        self.root.mainloop()
        
    def load_image(self):
        # clean the canvas before loading
        img_addr = join('features', "{gene_name}.png".format(gene_name = self.current_gene))
        self.img_data = cv2.imread(img_addr)
        self.img_data = cv2.cvtColor(self.img_data, cv2.COLOR_BGR2RGB)
        
        self.bg_img = Image.fromarray(self.img_data)
        self.photo = ImageTk.PhotoImage(image=self.bg_img)
        self.canvas.create_image(0, 0, image=self.photo, anchor=tk.NW)
        
        self.gene_name_label['text'] = self.current_gene
        
        self.write_gene_description()
        self.write_group_description()
        
        self.query_enter.delete(0, tk.END)
        
    def previous_group(self, event):
        self.current_cluster_index -= 1
        if self.current_cluster_index < 0:
            self.current_cluster_index = len(self.clusters) - 1
        self.current_cluster = self.clusters[self.current_cluster_index]
        self.current_gene_index = 0 
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.update_genes()
        self.group_list.selection_clear(0, tk.END)
        self.group_list.select_set(first = self.current_cluster_index)
        self.group_list.see(self.current_cluster_index)
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.gene_list.delete(0, tk.END)
        self.update_genes()
        self.load_image()
            
    def next_group(self, event):
        self.current_cluster_index += 1
        if self.current_cluster_index >= len(self.clusters):
            self.current_cluster_index = 0
        self.current_cluster = self.clusters[self.current_cluster_index]
        self.current_gene_index = 0 
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.update_genes()
        self.group_list.selection_clear(0, tk.END)
        self.group_list.select_set(first = self.current_cluster_index)
        self.group_list.see(self.current_cluster_index)
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.update_genes()
        self.load_image()
         
    def go_to_previous(self, event):
        self.current_gene_index -= 1
        if self.current_gene_index < 0:
            self.current_gene_index = len(self.data["Cluster_{number}".format(number = self.current_cluster)][0]) - 1
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.gene_name_label['text'] = self.current_gene
        self.gene_list.selection_clear(0, tk.END)
        self.gene_list.select_set(first = self.current_gene_index)
        self.gene_list.see(self.current_gene_index)
        self.load_image()
        
    def go_to_next(self, event):
        self.current_gene_index += 1
        if self.current_gene_index == len(self.data["Cluster_{number}".format(number = self.current_cluster)][0]):
            self.current_gene_index = 0
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.gene_name_label['text'] = self.current_gene
        self.gene_list.selection_clear(0, tk.END)
        self.gene_list.select_set(first = self.current_gene_index)
        self.gene_list.see(self.current_gene_index)
        self.load_image()
        
    def update_groups(self):
        for cluster in self.clusters:
            self.group_list.insert(tk.END, "Group_{number}".format(number = cluster))
        self.group_list.select_set(0)
    
    def update_genes(self):
        self.gene_list.delete(0, tk.END)
        for gene_name in self.data["Cluster_{number}".format(number = self.current_cluster)][0]:
            self.gene_list.insert(tk.END, gene_name)
        self.gene_list.select_set(0)
        self.gene_label['text'] = "Gene names ({numbers})".format(numbers = len(self.data["Cluster_{number}".format(number = self.current_cluster)][0]))
        
    def select_gene(self, event):
        self.current_gene_index = self.gene_list.curselection()[0]
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.load_image()
        
    def select_group(self, event):
        self.current_cluster_index = self.group_list.curselection()[0]
        self.current_cluster = self.clusters[self.current_cluster_index]
        self.current_gene_index = 0 
        self.current_gene = self.data["Cluster_{number}".format(number = self.current_cluster)][0][self.current_gene_index]
        self.gene_list.delete(0, tk.END)
        self.update_genes()
        self.load_image()
        
    def write_gene_description(self):
        gene_field      = self.gene_info.GeneSymbol == self.current_gene
        gene_symbol     = "Gene symbol: {sym}\n".format(sym = self.current_gene)
        if len(np.unique(gene_field)) == 1:
            gene_name       = "Gene name: Not available\n"
            gene_family     = "Gene family: Not available\n"
            reactom_pathway = "Reactom pathway: Not available\n"
            gene_summary    = "Gene summary: Not available\n"
        else:
            gene_name       = "Gene name: {x_factor}\n".format(x_factor = self.gene_info.GeneName[gene_field].values[0])
            gene_family     = "Gene family: {x_factor}\n".format(x_factor = self.gene_info.GeneFamily[gene_field].values[0])
            reactom_pathway = "Reactom pathway: {x_factor}\n".format(x_factor = self.gene_info.ReactomPathway[gene_field].values[0])
            gene_summary    = "Gene summary: {x_factor}\n".format(x_factor = self.gene_info.GeneSummary[gene_field].values[0])
        gene_field = "\n".join([gene_symbol, gene_name, gene_family, reactom_pathway, gene_summary])
        self.gene_description_content.delete('1.0', tk.END)
        self.gene_description_content.insert(tk.END, gene_field)
        
    def write_group_description(self):
        text_info = self.data["Cluster_{number}".format(number = self.current_cluster)][1]
        self.group_description_content.delete('1.0', tk.END)
        self.group_description_content.insert('1.0', text_info)
        
    def query_gene_name(self, event):
        query_entry = self.query_enter.get()
        if query_entry in self.clustering.GeneNames.values:
            self.current_gene = query_entry
            self.current_cluster = self.clustering.Cluster[self.clustering.GeneNames == self.current_gene].values[0]
            self.current_cluster_index = np.where(self.clusters == self.current_cluster)[0][0]
            self.current_gene_index = self.data["Cluster_{number}".format(number = self.current_cluster_index)][0]
            self.current_gene_index = self.current_gene_index.index(self.current_gene)
            self.update_genes()
            self.gene_list.selection_clear(0, tk.END)
            self.gene_list.select_set(first = self.current_gene_index)
            self.gene_list.see(self.current_gene_index)
            self.group_list.selection_clear(0, tk.END)
            self.group_list.select_set(first = self.current_cluster_index)
            self.group_list.see(self.current_cluster_index)
            self.load_image()
        else:
           self.query_enter.delete(0, tk.END)
           self.query_enter.insert(0, "Gene not found: {entrquery}".format(entrquery = query_entry))
    
    def update_group_description(self, event):
        text_info = self.group_description_content.get('1.0', tk.END)
        self.data["Cluster_{number}".format(number = self.current_cluster)][1] = text_info
        
    def save_changes(self, event):
        rmtree("group_descriptions")
        mkdir("group_descriptions")
        for cluster in self.clusters:
            cluster_name = "Cluster_{number}".format(number = cluster)
            fobj = open(join('group_descriptions', cluster_name), 'w')
            text_info = self.data[cluster_name][1]
            fobj.writelines(text_info)
            fobj.close()
        self.clustering.to_csv(self.clustering_addr)
            
    def changeGroupAssignment(self, event):
        self.new_window = tk.Toplevel(self.root)
        
        self.new_window.grab_set()
        
        self.options_frame = tk.Frame(self.new_window)
        self.actions_frame = tk.Frame(self.new_window)
        
        self.options = tk.Listbox(self.options_frame, exportselection = 0)
        self.options.config(selectmode = tk.SINGLE)
        self.options.pack()
        
        self.change_assignment_button  = tk.Button(self.actions_frame, text = 'Change assignment')
        self.show_selected_description = tk.Text(self.actions_frame, borderwidth=1, relief="solid", width = 28, height = 10)
        
        self.change_assignment_button.grid( row = 0, column = 0, sticky = tk.N)
        self.show_selected_description.grid(row = 2, column = 0, sticky = tk.N)
       
        for cluster in self.clusters:
            cluster_name = "Cluster_{number}".format(number = cluster)
            self.options.insert(tk.END, cluster_name)
        self.options.select_set(first = 0)
        self.options.see(0)
        
        self.show_selected_description.delete('1.0', tk.END)
        self.show_selected_description.insert('1.0', self.data['Cluster_0'][1])
        
        self.options_frame.grid(row = 0, column = 0, sticky = tk.N)
        self.actions_frame.grid(row = 0, column = 1, sticky = tk.N)
        
        self.new_window.bind("<Up>", self.assignment_up)
        self.new_window.bind("<Down>", self.assignment_down)
        self.change_assignment_button.bind("<ButtonRelease-1>", self.change_assignment)
        
    def assignment_up(self, event):
        asgn_index = self.options.curselection()[0]
        asgn_index -= 1
        if asgn_index < 0:
            asgn_index = len(self.clusters) - 1
        self.options.selection_clear(0, tk.END)
        self.options.select_set(first = asgn_index)
        self.options.see(asgn_index) 
        self.show_selected_description.delete('1.0', tk.END)
        self.show_selected_description.insert('1.0', self.data['Cluster_{number}'.format(number = asgn_index)][1])
    
    def assignment_down(self, event):
        asgn_index = self.options.curselection()[0]
        asgn_index += 1
        if asgn_index >= len(self.clusters):
            asgn_index =  0
        self.options.selection_clear(0, tk.END)
        self.options.select_set(first = asgn_index)
        self.options.see(asgn_index)
        self.show_selected_description.delete('1.0', tk.END)
        self.show_selected_description.insert('1.0', self.data['Cluster_{number}'.format(number = asgn_index)][1])
        
    def change_assignment(self, event):
        asgn_index = self.options.curselection()[0]
        self.clustering.Cluster[self.clustering.GeneNames == self.current_gene] = asgn_index
        self.data = {}
        for cluster in self.clusters:
            gene_names = [gene_name for gene_name in self.clustering.GeneNames[self.clustering.Cluster == cluster]]
            cluster_name = "Cluster_{number}".format(number = cluster)
            fobj = open(join('group_descriptions', cluster_name), 'r')
            cluster_description = fobj.read().strip()
            fobj.close()
            self.data[cluster_name] = [gene_names, cluster_description]
        self.query_enter.delete(0, tk.END)
        self.query_enter.insert(0, self.current_gene)
        self.query_gene_name(event)
        self.new_window.destroy()
        self.new_window.grab_release()
        
    def merge_groups(self, event):
        self.new_window = tk.Toplevel(self.root)
        
        self.new_window.grab_set()
        
        self.up_frame   = tk.Frame(self.new_window)
        self.down_frame = tk.Frame(self.new_window)
        
        self.up_frame.grid(  row = 0, column = 0)
        self.down_frame.grid(row = 1, column = 0)
        
        self.confirm_merging_button = tk.Button(self.up_frame, text = 'Confirm merging')
        self.confirm_merging_button.grid(row = 0, column = 0)
        
        self.set1_listbox     = tk.Listbox(self.down_frame, exportselection = 0)
        self.set1_description = tk.Text(self.down_frame, width = 28, height = 10)
        self.set2_listbox     = tk.Listbox(self.down_frame, exportselection = 0)
        self.set2_description = tk.Text(self.down_frame, width = 28, height = 10)
        
        self.set1_listbox.config(selectmode = tk.SINGLE)
        self.set2_listbox.config(selectmode = tk.SINGLE)
        
        self.set1_listbox.grid(    row = 0, column = 0)
        self.set1_description.grid(row = 0, column = 1)
        self.set2_listbox.grid(    row = 0, column = 2)
        self.set2_description.grid(row = 0, column = 3)
        
        for cluster in self.clusters:
            cluster_name = "Cluster_{number}".format(number = cluster)
            self.set1_listbox.insert(tk.END, cluster_name)
            self.set2_listbox.insert(tk.END, cluster_name)
            
        self.set1_listbox.select_set(first = self.current_cluster_index)
        self.set1_listbox.see(self.current_cluster_index)
        
        self.set2_listbox.select_set(first = 0)
        self.set2_listbox.see(0)
        
        self.set1_description.delete('1.0', tk.END)
        self.set1_description.insert('1.0', self.data['Cluster_{number}'.format(number = self.current_cluster_index)][1])
        
        self.set2_description.delete('1.0', tk.END)
        self.set2_description.insert('1.0', self.data['Cluster_0'][1])
        
        self.set1_listbox.bind(          "<ButtonRelease-1>", self.setListbox1)
        self.set2_listbox.bind(          "<ButtonRelease-1>", self.setListbox2)
        self.confirm_merging_button.bind("<ButtonRelease-1>", self.confirm_merging)
        
    def setListbox1(self, event):
        sel_idx = self.set1_listbox.curselection()[0]
        self.set1_description.delete('1.0', tk.END)
        self.set1_description.insert('1.0', self.data['Cluster_{number}'.format(number = sel_idx)][1])
        
        
    def setListbox2(self, event):
        sel_idx = self.set2_listbox.curselection()[0]
        self.set2_description.delete('1.0', tk.END)
        self.set2_description.insert('1.0', self.data['Cluster_{number}'.format(number = sel_idx)][1])
        
    def confirm_merging(self, event):
        member_1 = self.set1_listbox.curselection()[0]
        member_2 = self.set2_listbox.curselection()[0]
        to_erase = max(member_1, member_2)
        to_keep =  min(member_1, member_2)
        cluster1 = "Cluster_{number}".format(number = to_keep)
        cluster2 = "Cluster_{number}".format(number = to_erase)
        self.data[cluster1][0].extend(self.data[cluster2][0])
        self.data[cluster1][1] = "{a}\n{b}".format(a = self.data[cluster1][1], b = self.data[cluster2][1])
        new_data = {}
        for cluster in self.data.keys():
            cluster_idx = int(cluster.split('_')[1])
            if cluster_idx < to_erase:
                new_data[cluster] = self.data[cluster]
            elif cluster_idx == to_erase:
                continue
            else:
                new_cluster = "Cluster_{number}".format(number = cluster_idx - 1)
                new_data[new_cluster] = self.data[cluster]
        self.data = new_data
        ClusterCol, GeneCol = [], []
        for cluster in self.data.keys():
            cluster_idx = cluster.split("_")[1]
            genes = self.data[cluster][0]
            ClusterCol.extend(len(genes) * [cluster_idx, ])
            GeneCol.extend(genes)
        self.clustering = pd.DataFrame.from_dict({"Cluster": ClusterCol, "GeneNames": GeneCol})
        self.clustering['Cluster'] = pd.to_numeric(self.clustering['Cluster'])
        
        self.clusters = np.array([int(a) for a in self.clustering.Cluster.unique()])
        np.ndarray.sort(self.clusters)
        
        self.group_list.delete(0, tk.END)
        for cluster in self.clusters:
            self.group_list.insert(tk.END, "Group_{number}".format(number = cluster))
        self.group_list.select_set(first = to_keep)
        self.group_list.see(to_keep)
        self.current_cluster_index = int(to_keep)
        
        self.query_enter.delete(0, tk.END)
        self.query_enter.insert(0, self.current_gene)
        self.query_gene_name(event)
        
        self.new_window.destroy()
        self.new_window.grab_release()
        
feature_addrs = [join("features", addr) for addr in listdir('features') if addr[-3:] == "png"]
clus_addr      = "clustering.csv"
gene_info_addr = "gene_info.csv"

GeneViewer(feature_addrs, clus_addr, gene_info_addr)


    
    
