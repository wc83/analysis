#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 10:08:56 2018

@author: william
"""

# correlation for when polarisation is important

def corel_pi(wav1,wav2,shift):
    from obspy.signal.cross_correlation import correlate

    corell=correlate(wav1, wav2, shift, demean=True,normalize='naive',domain='time' )
    top=corell.argmax()
    top_v = corell[top]
    top = (shift)-top


    return(top_v,top,corell)