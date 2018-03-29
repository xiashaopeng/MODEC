# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 18:55:56 2018

@author: xsp
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import re

class modec_xml_opfile:
    # 初始化类
    def __init__(self, filename):
        # 文件名信息存储
        self.xml_opfile = filename
        
        # 各物理量单位字符串
        self.time_unit = ''
        self.power_unit = ''
        self.flux_unit = ''
        self.dens_unit = ''
        self.activity_unit = ''
        self.ampc_unit = ''
        self.wmpc_unit = ''
        self.toxi_unit = ''
        self.absrate_unit = ''
        self.prodrate_unit = ''
        
        # 各物理量ndarray
        self.time_array = np.array([])
        self.power_array = np.array([])
        self.flux_array = np.array([])
        self.kinf_array = np.array([])
        
        # 中子消失率和中子产生率的存储数组
        self.nucl_absrate_array = np.array([])
        self.nucl_prodrate_array = np.array([])
 
        # 放射性活度存储数组
        self.nucl_activity_array = np.array([])
        self.nucl_activity_loop_array = np.array([])
        self.nucl_activity_stockage_array = np.array([])
        
        # AMPC存储数组
        self.nucl_ampc_array = np.array([])
        self.nucl_ampc_loop_array = np.array([])
        self.nucl_ampc_stockage_array = np.array([])
        
        # 核素物质的量存储数组
        self.nucl_concentration_array = np.array([])
        self.nucl_concentration_loop_array = np.array([])
        self.nucl_concentration_stockage_array = np.array([])
        
        # 核素衰变热存储数组
        self.nucl_decayheat_array = np.array([])
        self.nucl_decayheat_loop_array = np.array([])
        self.nucl_decayheat_stockage_array = np.array([])

        # 核素放射性毒性存储数组
        self.nucl_toxi_array = np.array([])
        self.nucl_toxi_loop_array = np.array([])
        self.nucl_toxi_stockage_array = np.array([])
        
        # WMPC存储数组
        self.nucl_wmpc_array = np.array([])
        self.nucl_wmpc_loop_array = np.array([])
        self.nucl_wmpc_stockage_array = np.array([])
        
        # 核素名称，ID，原子序数以及相对原子质量
        self.nucl_name_array = np.array([])
        self.nucl_zai_array = np.array([],dtype=int)
        self.nucl_z_array = np.array([],dtype=int)
        self.nucl_a_array = np.array([])
        
	    # 解析XML文件
        self.parse_modec_xml_opfile()
		
    def parse_modec_xml_opfile(self):
        tree = ET.parse(self.xml_opfile)
        root = tree.getroot()
    
        for child in root:
            if child.tag == 'Time':
                self.time_unit = child.attrib["unit"]
                self.time_array = np.fromstring(child.text,sep=' ')
            if child.tag == 'Power':
                self.power_unit = child.attrib["unit"]
                self.power_array = np.fromstring(child.text,sep=' ')
            if child.tag == 'Flux':
                self.flux_unit = child.attrib["unit"]
                self.flux_array = np.fromstring(child.text,sep=' ')
            if child.tag == 'K-infinite':
                self.kinf_array = np.fromstring(child.text,sep=' ')
            if child.tag == 'Nuclides':
                zone_attrib = child.attrib["zone"]
                if zone_attrib == "core":
                    for childchild in child:
                        if childchild.tag == 'Concentration':
                            self.dens_unit = childchild.attrib["unit"]
                            # 核素浓度array
                            flag = 1
                            for childchildchild in childchild:
                                self.nucl_zai_array = np.append(self.nucl_zai_array,int(childchildchild.attrib["zai"]))
                                self.nucl_name_array = np.append(self.nucl_name_array,childchildchild.attrib["name"])
                                if flag == 1:
                                    self.nucl_concentration_array = np.append(self.nucl_concentration_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_concentration_array = np.vstack((self.nucl_concentration_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                
                        if childchild.tag == 'RadioActivity':
                            self.activity_unit = childchild.attrib["unit"]
                            # 放射性活度array
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_activity_array = np.append(self.nucl_activity_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_activity_array = np.vstack((self.nucl_activity_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                        
                        if childchild.tag == 'DecayHeat':
                            self.decayheat_unit = childchild.attrib["unit"]
                            # 衰变热
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_decayheat_array = np.append(self.nucl_decayheat_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_decayheat_array = np.vstack((self.nucl_decayheat_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'AMPC':
                            self.ampc_unit = childchild.attrib["unit"]
                            # AMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_ampc_array = np.append(self.nucl_ampc_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_ampc_array = np.vstack((self.nucl_ampc_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'WMPC':
                            self.wmpc_unit = childchild.attrib["unit"]
                            # WMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_wmpc_array = np.append(self.nucl_wmpc_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_wmpc_array = np.vstack((self.nucl_wmpc_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'RadioToxicity':
                            self.toxi_unit = childchild.attrib["unit"]
                            # 毒性
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_toxi_array = np.append(self.nucl_toxi_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_toxi_array = np.vstack((self.nucl_toxi_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                
                        if childchild.tag == 'NeuProdRate':
                            self.prodrate_unit = childchild.attrib["unit"]
                            # 中子产生率
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_prodrate_array = np.append(self.nucl_prodrate_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_prodrate_array = np.vstack((self.nucl_prodrate_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'NeuAbsRate':
                            self.absrate_unit = childchild.attrib["unit"]
                            # 中子吸收率
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_absrate_array = np.append(self.nucl_absrate_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_absrate_array = np.vstack((self.nucl_absrate_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                
                if zone_attrib == "stockage":
                    for childchild in child:
                        if childchild.tag == 'Concentration':
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_concentration_stockage_array = np.append(self.nucl_concentration_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_concentration_stockage_array = np.vstack((self.nucl_concentration_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                
                        if childchild.tag == 'RadioActivity':
                            # 放射性活度array
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_activity_stockage_array = np.append(self.nucl_activity_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_activity_stockage_array = np.vstack((self.nucl_activity_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                        
                        if childchild.tag == 'DecayHeat':
                            # 衰变热
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_decayheat_stockage_array = np.append(self.nucl_decayheat_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_decayheat_stockage_array = np.vstack((self.nucl_decayheat_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'AMPC':
                            # AMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_ampc_stockage_array = np.append(self.nucl_ampc_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_ampc_stockage_array = np.vstack((self.nucl_ampc_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'WMPC':
                            # WMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_wmpc_stockage_array = np.append(self.nucl_wmpc_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_wmpc_stockage_array = np.vstack((self.nucl_wmpc_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'RadioToxicity':
                            # 毒性
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_toxi_stockage_array = np.append(self.nucl_toxi_stockage_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_toxi_stockage_array = np.vstack((self.nucl_toxi_stockage_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                            
                if zone_attrib == "loop":
                    for childchild in child:
                        if childchild.tag == 'Concentration':
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_concentration_loop_array = np.append(self.nucl_concentration_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_concentration_loop_array = np.vstack((self.nucl_concentration_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                
                        if childchild.tag == 'RadioActivity':
                            # 放射性活度array
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_activity_loop_array = np.append(self.nucl_activity_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_activity_loop_array = np.vstack((self.nucl_activity_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                        
                        if childchild.tag == 'DecayHeat':
                            # 衰变热
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_decayheat_loop_array = np.append(self.nucl_decayheat_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_decayheat_loop_array = np.vstack((self.nucl_decayheat_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'AMPC':
                            # AMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_ampc_loop_array = np.append(self.nucl_ampc_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_ampc_loop_array = np.vstack((self.nucl_ampc_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'WMPC':
                            # WMPC
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_wmpc_loop_array = np.append(self.nucl_wmpc_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_wmpc_loop_array = np.vstack((self.nucl_wmpc_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1
                                            
                        if childchild.tag == 'RadioToxicity':
                            # 毒性
                            flag = 1
                            for childchildchild in childchild:
                                if flag == 1:
                                    self.nucl_toxi_loop_array = np.append(self.nucl_toxi_loop_array,np.fromstring(childchildchild.text,sep=' '))
                                else:
                                    self.nucl_toxi_loop_array = np.vstack((self.nucl_toxi_loop_array,np.fromstring(childchildchild.text,sep=' ')))
                                flag += 1   
        self.nucl_z_array = np.int_(self.nucl_zai_array / 10000)
        self.nucl_a_array =np.float_(np.int_(self.nucl_zai_array / 10) - self.nucl_z_array * 1000)

    def plt_time_evolution(self, nucl_name = [], **kwargs):
        # key = y_axis, x_axis 
        # value = 'linear', 'log'
        # key = format
        # value = 'sep' 'comb'
        
        # plt_array = []
        if 'object' in kwargs:
            value = kwargs['object']
            if value == 'concentration':
               plt_array = self.nucl_concentration_array
            elif value == 'activity':
                plt_array = self.nucl_activity_array
            elif value == 'decayheat':
                plt_array = self.nucl_decayheat_array
            elif value == 'toxicity':
                plt_array = self.nucl_toxi_array
            else:
                raise Exception('Unknown \'object\' value')
        else:
            plt_array = self.nucl_concentration_array

        if 'format' in kwargs:
            value = kwargs['format']
            if value == 'sep':
                for nucl in nucl_name:
                    if not any(char.isdigit() for char in nucl):
                        elem = [i for i, item in enumerate(self.nucl_name_array) if nucl == re.match(r'[a-zA-Z]+', item)[0]]
                        elem_concentration_array = np.full_like(plt_array[0], 0.0)
                        for elem_ in elem:
                            elem_concentration_array += plt_array[elem_]
                        plt.plot(self.time_array/3600./24.,elem_concentration_array,'o-')
                        plt.text(self.time_array.max() / 3600. /24. ,
                        elem_concentration_array.max() * .8 + elem_concentration_array.min() * .2, nucl)
                    else:
                        index = np.where(self.nucl_name_array == nucl)
                        plt.plot(self.time_array / 3600. / 24., plt_array[index[0][0]],'o-')
                        plt.text(self.time_array.max() / 3600. / 24. * 0.8,
                        plt_array[index[0][0]].max() * .8 + 
                        plt_array[index[0][0]].min() * .2
                        , nucl)
                    
                    if 'x_axis' in kwargs:
                        plt.xscale(kwargs['x_axis'])
                    if 'y_axis' in kwargs:
                        plt.yscale(kwargs['y_axis'])

                    plt.ylabel(kwargs['object'])
                    plt.xlabel("time(days)")
                    plt.grid(True)
                    # plt.show()
                    plt.savefig(nucl)
                    plt.close()

            elif value == 'comb':
                com_nucl = np.full_like(plt_array[0], 0.0)
                for nucl in nucl_name:
                    if not any(char.isdigit() for char in nucl):
                        elem = [i for i, item in enumerate(self.nucl_name_array) if nucl == re.match(r'[a-zA-Z]+', item)[0]]
                        for elem_ in elem:
                            com_nucl += plt_array[elem_]
                    else :
                        index = np.where(self.nucl_name_array == nucl) 
                        com_nucl += plt_array[index[0][0]]
                plt.plot(self.time_array / 3600. / 24., com_nucl, 'o-')
                plt.text(self.time_array.max() / 3600. / 24. * 0.8,
                        com_nucl.max() * .8 + com_nucl.min() * .2
                        , nucl)

                if 'x_axis' in kwargs:
                    plt.xscale(kwargs['x_axis'])
                if 'y_axis' in kwargs:
                    plt.yscale(kwargs['y_axis'])

                plt.ylabel(kwargs['object'])
                plt.xlabel("time(days)")
                plt.grid(True)
                # plt.show()
                if 'savefig' in kwargs:
                    plt.savefig(kwargs['savefig'])
                else:
                    plt.savefig('0')
                # plt.close()

        else:
            plt.figure()
            for nucl in nucl_name:
                if not any(char.isdigit() for char in nucl):
                    elem = [i for i, item in enumerate(self.nucl_name_array) if nucl == re.match(r'[a-zA-Z]+', item)[0]]
                    elem_concentration_array = np.full_like(plt_array[0], 0.0)
                    for elem_ in elem:
                        elem_concentration_array += plt_array[elem_]
                    plt.plot(self.time_array/3600./24.,elem_concentration_array,'o-')
                    plt.text(self.time_array.max() / 3600. /24. ,
                    elem_concentration_array.max() * .8 + elem_concentration_array.min() * .2, nucl)
                else :
                    index = np.where(self.nucl_name_array==nucl)
                    plt.plot(self.time_array/3600/24,plt_array[index[0][0]],'o-')
                    plt.text(self.time_array.max()/3600/24*0.8, 
                        plt_array[index[0][0]].max()*.8 + 
                        plt_array[index[0][0]].min()*.2, nucl)
                #                ,bbox=dict(facecolor='red', alpha=0.5))
        #            plt.annotate(nucl,
        #                         xy = (self.time_array.max()/3600/24*0.8, 
        #                 self.nucl_concentration_array[index[0][0]].max()*.8 + 
        #                self.nucl_concentration_array[index[0][0]].min()*.2),
        #                xytext = (self.time_array.max()/3600/24*0.8, 
        #                 self.nucl_concentration_array[index[0][0]].max()*.8 + 
        #                self.nucl_concentration_array[index[0][0]].min()*.2))
        
            if 'x_axis' in kwargs:
                plt.xscale(kwargs['x_axis'])
            if 'y_axis' in kwargs:
                plt.yscale(kwargs['y_axis'])

            plt.ylabel(kwargs['object'])
            plt.xlabel("time(days)")

            plt.grid(True)
            if 'savefig' in kwargs:
                plt.savefig(kwargs['savefig'])
            else:
                plt.savefig('0')
            # plt.show()

    def output_time_evolution(self, nucl_name = [], **kwargs):
        # key = y_axis, x_axis 
        # value = 'linear', 'log'
        # key = format
        # value = 'sep' 'comb'
        
        # plt_array = []
        if 'object' in kwargs:
            value_o = kwargs['object']
            if value_o == 'concentration':
                nucl_array = self.nucl_concentration_array
                nucl_str = self.dens_unit
            elif value_o == 'activity':
                nucl_array = self.nucl_activity_array
                unit_str = self.activity_unit
            elif value_o == 'decayheat':
                nucl_array = self.nucl_decayheat_array
                unit_str = self.decayheat_unit
            elif value_o == 'toxicity':
                nucl_array = self.nucl_toxi_array
                unit_str = self.toxi_unit
            else:
                raise Exception('Unknown \'object\' value')
        else:
            value_o = 'concentration'
            nucl_array = self.nucl_concentration_array
            nucl_str = self.dens_unit

        time_ = self.time_array/3600./24./365.25
            
        for nucl in nucl_name:
            # 文件名
            if nucl == "Total":
                # 文件名
                out_file = open('Total_' + value_o + "_evolution.dat",'w')
                
                # 文件头
                out_file.write("========> {object_} evolution of total nuclides ({unit_}) <========\n\n".
                               format(object_ = value_o, unit_ = unit_str))
                # 输出结果
                out_file.write("{Time:<15s}".format(Time= 'Time(years)'))
    
                elem_concentration_array = np.full_like(nucl_array[0], 0.0)
                
                for i in range(nucl_array.shape[0]):
                    elem_concentration_array += nucl_array[i]
                    
                out_file.write("{nuclide_name:<16s}".format(nuclide_name = nucl))
                out_file.write('\n')
                
                for i in range(time_.size):
                    out_file.write("{:<15E}".format(time_[i]))
                    out_file.write("{:<16.9E}".format(elem_concentration_array[i]))
                    out_file.write('\n')
                
                out_file.close()
            else:
                out_file = open(nucl + '_' + value_o + "_evolution.dat",'w')
                
                # 文件头
                out_file.write("========> {object_} evolution of {nucl_}({unit_}) <========\n\n".
                               format(object_ = value_o, nucl_ = nucl, unit_ = unit_str))
                # 输出结果
                out_file.write("{Time:<15s}".format(Time= 'Time(years)'))
                if not any(char.isdigit() for char in nucl):
                    elem = [i for i, item in enumerate(self.nucl_name_array) if nucl == re.match(r'[a-zA-Z]+', item)[0]]
                    elem_concentration_array = np.full_like(nucl_array[0], 0.0)
                
                    for elem_ in elem:
                        out_file.write("{nuclide_name:<16s}".format(nuclide_name = self.nucl_name_array[elem_]))
                        elem_concentration_array += nucl_array[elem_]
                    out_file.write("{nuclide_name:<16s}".format(nuclide_name="Total"))
                    out_file.write('\n')
                    
                    for i in range(time_.size):
                        out_file.write("{:<15E}".format(time_[i]))
                        for elem_ in elem:
                            out_file.write("{:<16.9E}".format(nucl_array[elem_][i]))
                        out_file.write("{:<16.9E}".format(elem_concentration_array[i]))
                        out_file.write('\n')
                else:
                    index = np.where(self.nucl_name_array == nucl)
                    out_file.write("{nuclide_name:<16s}".format(nuclide_name = self.nucl_name_array[index[0][0]]))
                    out_file.write('\n')
                    for i in range(time_.size):
                        out_file.write("{:<15E}".format(time_[i]))
                        out_file.write("{:<16.9E}".format(self.nucl_name_array[index[0][0]][i]))
                        out_file.write('\n')
                
                out_file.close()
