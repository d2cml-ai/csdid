# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 16:17:07 2023

@author: Carlos Guevara
"""

import matplotlib.pyplot as plt

# get_ipython().run_line_magic('matplotlib', 'qt') # To aopen separate window
# get_ipython().run_line_magic('matplotlib', 'inline') # In line graph

def gplot(ssresults, ax, ylim=None, xlab=None, ylab=None, title="Group", xgap=1,
           legend=True, ref_line=0, theming=True):
    if ylab is None:
        ylab = 'ATT'
    
    ssresults = ssresults[ssresults['year'].notnull()].copy()
    ssresults.loc[:, 'year'] = ssresults['year'].astype(int).astype(str)
    
    pre_points = ssresults.loc[ssresults['post'] == 0]
    post_points = ssresults.loc[ssresults['post'] == 1]
    
    ax.errorbar(pre_points['year'], pre_points['att'], yerr=pre_points['c']*pre_points['att_se'],
                 fmt='o', markersize=5, color='#e87d72', ecolor='#e87d72', capsize=5, label='Pre')   
    
    ax.errorbar(post_points['year'], post_points['att'], yerr=post_points['c']*post_points['att_se'],
                 fmt='o', markersize=5, color='#56bcc2', ecolor='#56bcc2', capsize=5, label='Post')  
    
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.set_xlabel(xlab)    
    ax.set_ylabel(ylab)    

    handles, labels = ax.get_legend_handles_labels()    
    
    if ref_line is not None:
        ax.axhline(ref_line, linestyle='dashed', color='#1F1F1F')
    if theming:
        ax.set_facecolor('white')
        ax.set_title(title, color="#1F1F1F", fontweight="bold", fontsize=10)
        ax.spines['bottom'].set_color('#1F1F1F')
        ax.spines['left'].set_color('#1F1F1F')
        ax.tick_params(axis='x', colors='#1F1F1F')
        ax.tick_params(axis='y', colors='#1F1F1F')
        if not pre_points.empty and not post_points.empty:
            ax.legend(handles[0:2], labels[0:2], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))
        elif not pre_points.empty:
            ax.legend(handles[:1], labels[:1], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))
        elif not post_points.empty:
            ax.legend(handles[1:2], labels[1:2], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))     
    if not legend:
        ax.legend().set_visible(False)
        
    return ax


def splot(ssresults, ax, ylim=None, xlab=None, ylab=None, title="Group",
          legend=True, ref_line=0, theming=True):
    
    if xlab is None:
        xlab = 'Group'
    if ylab is None:
        ylab = 'ATT'

    ssresults['year'] = ssresults['year'].copy().astype(str)
    
    pre_points = ssresults.loc[ssresults['post'] == 0]
    post_points = ssresults.loc[ssresults['post'] == 1]
    
    ax.errorbar(pre_points['year'], pre_points['att'], yerr=pre_points['c']*pre_points['att_se'],
                 fmt='o', markersize=5, color='#e87d72', ecolor='#e87d72', capsize=5, label='Pre')   
    
    ax.errorbar(post_points['year'], post_points['att'], yerr=post_points['c']*post_points['att_se'],
                 fmt='o', markersize=5, color='#56bcc2', ecolor='#56bcc2', capsize=5, label='Post') 
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)

    handles, labels = ax.get_legend_handles_labels()    
    
    if ylim is not None:
        ax.set_ylim(ylim)
    
    if ref_line is not None:
        ax.axhline(ref_line, linestyle='dashed', color='#1F1F1F')
    
    if theming:
        ax.set_facecolor('white')
        ax.set_title(title, color="#1F1F1F", fontweight="bold", fontsize=12)
        ax.spines['bottom'].set_color('#1F1F1F')
        ax.spines['left'].set_color('#1F1F1F')
        ax.tick_params(axis='x', colors='#1F1F1F')
        ax.tick_params(axis='y', colors='#1F1F1F')
        if not pre_points.empty and not post_points.empty:
            ax.legend(handles[0:2], labels[0:2], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))
        elif not pre_points.empty:
            ax.legend(handles[:1], labels[:1], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))
        elif not post_points.empty:
            ax.legend(handles[1:2], labels[1:2], loc='lower center',fontsize='small', ncol=2, bbox_to_anchor=(0.5,-0.27))   
            
    if not legend:
        ax.legend().set_visible(False)
    
    return ax