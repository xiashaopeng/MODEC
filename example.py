# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 19:23:17 2018

@author: xsp
"""

from modec_xml_read import modec_xml_opfile

xml_file = modec_xml_opfile("input_out.xml")
#xml_file.parse_modec_xml_opfile()
#print(xml_file.time_array.shape)
print(xml_file.nucl_zai_array[100])
print(xml_file.nucl_z_array[100])
print((xml_file.nucl_a_array[100]))

# 画图示例
xml_file.plt_time_evolution(("Th232","U233","Xe"), y_axis = 'log', x_axis = 'linear', format = 'sep')
