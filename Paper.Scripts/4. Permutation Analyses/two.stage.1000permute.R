#!/bin/bash

# 11.3.2020
# This script will generate 1000 permutation orders for all tissues

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
        stop("One argument must be supplied (fileID).", call.=FALSE)
}

# load libs
library(analyzeR)
library(yangR)
library(magrittr)

# get the tissue
num<-args[1]
num <- as.numeric(num)
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')
char <- all_tissues[num]

setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# load dataframe

if (!file.exists('full.frame.rds')){print(paste0('full.frame.rds does not exist for ', char))} else{
        load('full.frame.rds')

	order.for.permutes <- as.data.frame(matrix(nrow = nrow(full.frame.ct.corr), ncol = 1000))
	set.seed(1)
	for(i in 1:1000)
	{
		order.for.permutes[,i] <- sample(full.frame.ct.corr$SUBJID)
	}

	save(order.for.permutes, file =  paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/order.for.permutes1000.rds'))
}

