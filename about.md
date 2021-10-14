---
layout: page
title: About
permalink: /about/
---

# Welcome to the Umit (*U*niversal *m*esoscale *I*maging *T*oolbox) Wiki!  
## About  
Recent advances in genetic tools, such as fast-responding calcium indicators (e.g., GCaMP6f) as well as the increasingly accessible high-speed, high-resolution and high throughput longitudinal recording systems allowed the acquisition of faster and more detailed images at mesoscale level for months on large cohorts of animals. This generates large amounts of data that need to be stored, pre-processed and analyzed in a timely manner which can become a limiting factor as experiments become larger and more complex.   

The *umIToolbox*, or *umIT* for short,  was created to provide a set of tools to manage and automate critical steps of processing large-scale imaging datasets. The toolbox is currently in active development and new features are added periodically. *umIT* is written in MATLAB and, currently, it consists mainly of four modules that manages the experiment metadata (i.e., information about experimental subjects and recordings), creates pipelines to import and pre-process data, provides data visualization tools and basic group comparison and statistical analysis (*currently under construction). All modules and functionalities can be accessed through the *umIToolbox main GUI*.   
[[https://github.com/S-Belanger/Umit/blob/master/wiki_images/umIT_4axis.png|alt=umIT_concept1]]
###### Fig 1. The 4 main roles of **umIT**.

***



## How the toolbox manages data   
A typical imaging **project** consists of one or more cohorts of **subject**s (e.g. mice) that undergo one or more **acquisition** (i.e. recording) sessions. Frequently, other recording **modalities** are associated with the imaging data such as behavioral responses, eye/body tracking, etc. The toolbox follows the same organization principle where one can manage subjects, acquisitions and recording modalities for a given project.   

The management of the project datasets is done through a set of objects* that control the project (class Protocol), subjects (class Subject), acquisitions (class Acquisition) and modalities (abstract class modality and sub classes such as *FluorescenceImaging* class, for imaging data). Once the project file is created, the toolbox can automatically detect new subjects or acquisitions and update the file.   

In addition to the experiment management tool, the toolbox provides a series of **analysis functions** and a **pipeline** engine that allows the automation of several steps of the data processing.    

These two modules (experiment management and data analysis) were created to be adaptable and easily extensible in order to fulfill the needs of different neuroimaging experiment designs.    
   

[[https://github.com/S-Belanger/Umit/blob/master/wiki_images/umIT_concept_org_img_exp.png|alt=umIT_concept2]]   
###### Fig 2. Conceptual organization of a typical imaging project.


