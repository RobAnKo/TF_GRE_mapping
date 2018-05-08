#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 18:46:05 2018

@author: rkoch
"""
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO



records = list(SeqIO.parse('/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters.fna', 'fasta'))

record = records[0]

p1 = [100]
p2 = [x +24 for x in p1]


for p,pp, i in zip(p1,p2, range(1,len(p1)+1)):
    record.features.append(SeqFeature(location = FeatureLocation(p, pp, strand = +1), type = 'GRE', id = 'G{}'.format(i)))

p3 = [180]
p4 = [x + 15 for x in p3]

for p,pp, i in zip(p3,p4, range(1,len(p3)+1)):
    record.features.append(SeqFeature(location = FeatureLocation(p, pp, strand = -1), type = 'TF', id = 'T{}'.format(i)))


p5 = []
p6 = []
for p,pp in zip(p1,p2):
    for ppp,pppp in zip(p3,p4):
        if ppp < p < pppp:
            p5.append(p)
            p6.append(pppp)
        elif p < ppp < pp:
            p5.append(ppp)
            p6.append(pp)
                
            
for p,pp, i in zip(p5,p6,range(1,len(p5)+1)):
    record.features.append(SeqFeature(location = FeatureLocation(p, pp, strand = None), type = 'overlap', id = 'o{}'.format(i)))

gd_diagram = GenomeDiagram.Diagram("Promoter region with GRE and TF binding")
gd_track_for_features = gd_diagram.new_track(1, name="GREs and TFBS")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type == "GRE":
        color = colors.green
        gd_feature_set.add_feature(feature, color = color, label = feature.id, label_size=12, label_angle=20, sigil = "OCTO")
    if feature.type == "TF":
        color = colors.red
        gd_feature_set.add_feature(feature, color = color, label = feature.id, label_size=12, label_angle=20, sigil = "OCTO")
#    else:
#        gd_feature_set.add_feature(feature, color = colors.yellow, label = False)
    
gd_diagram.draw(format="linear", orientation="landscape", pagesize=(300,80),
                fragments=1, start=0, end=len(record))


gd_diagram.write("both2_pos.eps", "eps")