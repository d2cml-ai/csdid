# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 16:17:07 2023

@author: Carlos Guevara
"""

from plotnine import ggplot, aes, geom_point, geom_errorbar, scale_y_continuous, scale_x_continuous, \
    scale_color_manual, labs, geom_hline, theme, element_text, element_rect, element_line, \
    scale_x_discrete

# get_ipython().run_line_magic('matplotlib', 'qt') # To aopen separate window
# get_ipython().run_line_magic('matplotlib', 'inline') # In line graph

def gplot(ssresults, ylim=None, xlab=None, ylab=None, title="Group", xgap=1,
           legend=True, ref_line=0, theming=True):
    if ylab is None:
        ylab = 'ATT'
    dabreaks = ssresults['year'][::xgap]
    ssresults['post'] = ssresults['post'].astype('category')

    post_levels = list(ssresults['post'].unique().categories)
    
    if 0 in post_levels:
        color_values = ["#e87d72", "#56bcc2"]
        label_values = ['Pre', 'Post']
    else:
        color_values = ["#56bcc2"]
        label_values = ['Post']
    
    p = ggplot(ssresults,
               aes(x='year', y='att', ymin='att - c * att_se', ymax='att + c * att_se')) + \
        geom_point(aes(colour='post'), size=2) + \
        geom_errorbar(aes(colour='post'), width=0.15 , size=0.8 ) + \
        scale_y_continuous(limits=ylim) + \
        scale_x_continuous(breaks=list(dabreaks), labels=list(map(str, dabreaks))) + \
        scale_color_manual(drop=False, values=color_values, breaks=post_levels, labels=label_values) + \
        labs(x=xlab, y=ylab, title=title, color='') 
        
    if ref_line is not None:
        p += geom_hline(aes(yintercept=ref_line), linetype='dashed')

    if theming:
        p += theme(
                panel_background=element_rect(fill='white'),
                plot_title=element_text(color='#1F1F1F', fontweight='bold', size=10),
                axis_text=element_text(color='#1F1F1F'),
                axis_line=element_line(color='#1F1F1F'),
                strip_background=element_rect(fill='white'),
                legend_position=(.55, -.025),
                strip_text=element_text(color='#1F1F1F', fontweight='bold', size=9)
        ) 

    if not legend:
        p += theme(legend_position='none')

    return p


def splot(ssresults, ylim=None, xlab=None, ylab=None, title="Group",
          legend=True, ref_line=0, theming=True):
    if ylab is None:
        ylab = 'ATT'
    if xlab is None:
        xlab = 'Group'
    ssresults['year'] = ssresults['year'].astype('category')
    ssresults['post'] = ssresults['post'].astype('category')
    
    post_levels = list(ssresults['post'].unique().categories)  
   
    if 0 in post_levels:
        color_values = ["#e87d72", "#56bcc2"]
        label_values = ['Pre', 'Post']
    else:
        color_values = ["#56bcc2"]
        label_values = ['Post']

    p = ggplot(ssresults,
               aes(y='att', x='year', ymin='att - c * att_se', ymax='att + c * att_se')) + \
        geom_point(aes(colour='post'), size=2) + \
        geom_errorbar(aes(colour='post'), width=0.15 , size=0.8) + \
        scale_y_continuous(limits=ylim) + \
        scale_x_discrete(breaks=list(ssresults['year']))  + \
        scale_color_manual(drop=False, values=color_values,breaks=list(post_levels) , labels=label_values) + \
        labs(x=xlab, y=ylab, title=title, color='')
    
    if ref_line is not None:
        p += geom_hline(aes(yintercept=ref_line), linetype='dashed')
    
    if theming:
        p += theme(
                panel_background=element_rect(fill='white'),
                plot_title=element_text(color='#1F1F1F', fontweight='bold', size=10),
                axis_text=element_text(color='#1F1F1F'),
                axis_line=element_line(color='#1F1F1F'),
                strip_background=element_rect(fill='white'),
                legend_position=(.55, -.025)
            )
    
    if not legend:
        p += theme(legend_position='none')
    
    return p