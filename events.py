from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import medfilt
from scipy.stats import median_abs_deviation as mad
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy import signal
import bokeh.plotting as bpl
import pandas as pd
import numpy as np
import pickle

def EventForm(xdata, a, b, k, t0, ton, toff):
    return [k*(t - t0) + b if t < t0 else a*(1 - np.exp((t0 - t)/ton))*np.exp((t0 - t)/toff) + b for t in xdata]

def BckgLine(xdata, k, b):
    return [k*x + b for x in xdata]



def clnm(num):
    #Returns color for each number as in Moscow Metro
    return {
    1:"red",        
    2:"green",      
    3:"mediumblue",        
    4:"cyan",        
    5:"sienna", 
    6:"darkorange",    
    7:"mediumvioletred",      
    8:"gold",   
    9:"magenta",
    0:"lawngreen"}.get(num%10) 


def SearchEventsByLinearCross(time, traces, opts):
    #search for significant events by double thresholding (by peak amplitude and by the peak of derivative).
    #the precise time of the event is defined as crossing between the derivative line and the backgrounfd line 

    def BiLinearFit(xdata, k, t0): #a, fdp are to be nailed
        b  = fdp - a*(xdata[-1] - t0)
        return [k*(t - t0) + b if t < t0 else a*(t - t0) + b for t in xdata]

    
    events = []    
    if opts['draw']:
        bpl.output_file(opts['output_name'])
        p = bpl.figure(title = opts['output_name'], height = 1000, width = 1800)#, width_policy = 'fit')
    
    for cell_num, trace in enumerate(traces):
        #smoothing of the trace
        sm_trace = gaussian_filter1d(trace, sigma=opts['sigma'])
        #derivative calculation
        d_trace = np.gradient(sm_trace, 0.01)

        #calculation of thresholds relative to MAD
        thr_ampl = opts['thr_ampl']*mad(trace)
        thr_der = opts['thr_der']*mad(d_trace)  
        
        #peak searching
        d_peaks = signal.argrelmax(d_trace)[0]
        x_peaks = signal.argrelmax(sm_trace)[0]
                
        if opts['draw']:
            p.line(time, trace/np.max(trace) + cell_num, line_color = clnm(cell_num))
            p.line(time, d_trace/np.max(d_trace) + cell_num, line_color = clnm(cell_num), line_dash = 'dashed', line_alpha = 0.7)
            p.line(time, sm_trace/np.max(trace) + cell_num, line_color = clnm(cell_num), line_dash = 'dotted', line_alpha = 0.7)
            p.scatter(time[x_peaks], trace[x_peaks]/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), fill_alpha = 0.5, marker = 'inverted_triangle', size = 8)
            p.scatter(time[d_peaks], trace[d_peaks]/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), fill_alpha = 0.5, marker = 'diamond_dot', size = 8)
        
        evs = []  #list of this cell's events
        for idpk, d_peak in enumerate(d_peaks):
            x_left = max(0, d_peak - opts['fit_wnd'])
            #applying of the additional threshold by the peak amplitude
            x_peak = len(trace)-1 if all (x_peaks < d_peak) else x_peaks[x_peaks >= d_peak][0] #x_peak is the next peak to the d_peak 
            if idpk + 1 < len(d_peaks) and x_peak > d_peaks[idpk + 1] or x_peak == d_peak or d_peak <= 1 or d_peak >= len(trace)-3: #if there is no x_peak nearer than the next d_peak, skip this d_peak
                continue
            ampl = max(trace[d_peak:x_peak+1]) - min(trace[x_left:d_peak])
            if d_trace[d_peak] >= thr_der and ampl >= thr_ampl:
                #estimated values and bounds of fitting params; a, fdp, k, t0; a and fdp are nailed to f'(d_peak) and f(d_peak) respectively.
                #fitting function is BiLinearFit(xdata, a, fdp, k, t0)
                #local derivative calculation
                d_trace_precise = np.gradient(trace[d_peak-2:d_peak+2], time[d_peak-2:d_peak+2])
                a = max([max(d_trace_precise[1:-1]), d_trace[d_peak]])
                #a = (trace[d_peak+1] - trace[d_peak-1])/(time[d_peak+1] - time[d_peak-1])
                fdp = trace[d_peak]
                p0 = (0, (time[x_left] + time[d_peak])/2)
                bounds=((-np.inf, time[x_left]), (opts['max_k']*mad(trace)/opts['fit_wnd'], time[d_peak])) #((l,o,w,e,r,_), (h,i,g,h,e,r))
                try:
                    popt,_ = curve_fit(BiLinearFit, time[x_left:d_peak+1], trace[x_left:d_peak+1], p0 = p0, bounds = bounds) 
                except:
                    print(f'FAILED to detect event at cell {cell_num} time {time[d_peak]} s')
                    continue
                evs.append(dict(zip(['cell_num', 'd_peak','x_left', 'x_peak', 'ampl', 'a', 'fdp', 'k', 't0'],[cell_num, d_peak, x_left, x_peak, ampl, a, fdp, *popt])))
                print(f'Event detected: {evs[-1]}') 
                if opts['draw']:
                    fit = np.array(BiLinearFit(time[x_left:d_peak+1], *popt))
                    p.scatter(popt[1], -0.1 + cell_num, line_color = None, fill_color = clnm(cell_num), size = 5)
                    p.line(time[x_left:d_peak+1], fit/np.max(trace)+cell_num, line_color = clnm(cell_num), line_width =3.0, line_alpha = 0.5)
        events.append(evs) 
    if opts['draw']:
        bpl.show(p)
    return events
            
    


def SearchEventsDerivativeFit(time, traces, opts):
    #searches the events within a given short time window (~1s); possibly, it can be not enough peaks/pits there
    #additional opts: 
    #thr_der = threshold by derivative
    #fit_window_left = left window size for the performance of the fit, typically 2-3 frames
    #fit_window_right = right window size for the performance of the fit
    events = [] 
    if opts['draw']:
        bpl.output_file(opts['output_name'])
        p = bpl.figure(title = opts['output_name'], height = 1000, width = 1800)#, width_policy = 'fit')
    
    for cell_num, trace in enumerate(traces):
        #smoothing of the trace
        sm_trace = gaussian_filter1d(trace, sigma=opts['sigma'])
        #derivative calculation
        d_trace = np.gradient(sm_trace, 0.01)

        #calculation of thresholds relative to MAD
        thr = opts['thr_peak']*mad(trace)
        thr_der = opts['thr_der']*mad(d_trace)  
        
        #peak searching
        d_peaks = signal.argrelmax(d_trace)[0]
        x_peaks = signal.argrelmax(sm_trace)[0]
        x_pits = signal.argrelmin(sm_trace)[0]  

        if opts['draw']:
            p.line(time, trace/np.max(trace) + cell_num, line_color = clnm(cell_num))
            p.line(time, d_trace/np.max(d_trace) + cell_num, line_color = clnm(cell_num), line_dash = 'dashed', line_alpha = 0.7)
            p.line(time, sm_trace/np.max(trace) + cell_num, line_color = clnm(cell_num), line_dash = 'dotted', line_alpha = 0.7)
            p.scatter(time[x_peaks], trace[x_peaks]/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), fill_alpha = 0.5, marker = 'inverted_triangle', size = 8)
            p.scatter(time[d_peaks], trace[d_peaks]/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), fill_alpha = 0.5, marker = 'diamond_dot', size = 8)
            p.scatter(time[x_pits], trace[x_pits]/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), fill_alpha = 0.5, marker = 'triangle', size = 8)
            
        evs = []  #list of this cell's events
        for idpk, d_peak in enumerate(d_peaks):
            x_peak = len(trace)-1 if all (x_peaks < d_peak) else x_peaks[x_peaks > d_peak][0] #x_peak is the next peak to the d_peak 
            x_pit = len(trace)-1 if all (x_pits < x_peak) else x_pits[x_pits > x_peak][0]  #x_pit is the next pit to the x_peak 
            
            x_left = max(0, d_peak - opts['fit_window_left'])
            x_right = min(len(trace)-1, x_peak + opts['fit_window_right'], x_pit)
            if idpk + 1 < len(d_peaks): #this is to avoid double-counting in case of absence of x_peak between two d_peaks
                x_right = min(x_right, d_peaks[idpk + 1])
            if trace[x_peak] - np.min(trace[x_left:x_peak]) >= thr and d_trace[d_peak] >= thr_der:  #thresholding by peak and by derivative value simultaneously
                #estimated values and bounds of fitting params; a,b,k,t0,ton,toff
                p0 = (thr, trace[x_left], 0, time[d_peak], opts['est_ton'], opts['est_toff'])
                bounds=((thr, -np.inf, -np.inf, time[x_left], 0, 0), (np.inf, np.inf, opts['max_k']*thr/opts['est_ton'], time[x_peak], np.inf, np.inf)) #((l,o,w,e,r,_),(h,i,g,h,e,r))
                try:
                    popt,_ = curve_fit(EventForm, time[x_left:x_right], trace[x_left:x_right], p0 = p0, bounds = bounds)  
                except:
                    try:
                        p0 = (thr, trace[x_left], 0, time[x_left], opts['est_ton'], opts['est_toff'])
                        popt,_ = curve_fit(EventForm, time[x_left:x_right], trace[x_left:x_right], p0 = p0, bounds = bounds)
                    except:
                        print(f'FAILED to detect event at cell {cell_num} time {time[x_peak]} s')
                        continue
                fit = EventForm(time[x_left:x_right], *popt)
                ampl = np.max(fit) - popt[1]  #relative amplitude of the event
                evs.append(dict(zip(['cell_num','ampl','a','b','k','t0','ton','toff','x_left','x_right'],[cell_num,ampl,*popt,x_left,x_right])))
                print(f'Event detected: {evs[-1]}')
                if opts['draw']:
                    p.scatter(popt[3], np.max(fit)/np.max(trace) + cell_num, line_color = None, fill_color = clnm(cell_num), size = 5)
                    p.line(time[x_left:x_right], fit/np.max(trace) + cell_num, line_color = clnm(cell_num), line_width =3.0, line_alpha = 0.5)
        events.append(evs)  
    if opts['draw']:
        bpl.show(p)
    return events

      
def DrawSpEvents(tr_fname, sp_fname):
    #file reading
    bpl.output_notebook()
    traces = np.genfromtxt(tr_fname, delimiter = ',', skip_header = 1)[:,1:].T #note the transposition for better iterability
    sp_events = np.genfromtxt(sp_fname, delimiter = ',', skip_header = 1)[:,1:].T
    time = np.genfromtxt(tr_fname, delimiter = ',', skip_header = 1)[:,0]

    p = bpl.figure(title = tr_fname.split('\\')[-1], width = 1000)#, width_policy = 'fit')
    for cell_num, (trace, spikes) in enumerate(zip(traces, sp_events)):
        p.line(time, trace/np.max(trace) + cell_num, line_color = clnm(cell_num))
        p.scatter(time[spikes>0], cell_num - 0.1, line_color = None, fill_color = clnm(cell_num), size = 5)
    bpl.show(p)
