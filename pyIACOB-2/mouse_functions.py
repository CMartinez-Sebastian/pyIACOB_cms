import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/cms/Desktop/PhD/0_PROJECT/pyIACOB-2')
import spec as sp 
import time


##### Functions definition #####
coords,window,masks = [],[],[]

def onclick(event):
    '''
    GOAL : Function to mask interactively centain regions in the observed spectra

    Instructions
    ----------
        Interactive selection of the wavelengths to be masked. To select the limits
        to be included in coords, double-click on it with the left mouse bottom.
        When selecting an area, it will appear as grey-shadow band in the interactive
        spectrum plot.
        
        To correct a selected area, remove it by pressing right mouse bottom over
        the area to remove.

    Returns
    -------
    coords : List,
        List composed of several 2D lists for constructing the masks (See masks 
        function)

    '''
    
    
    
    if event.button == 1 and event.dblclick == True:
        global ix, iy
        ix = event.xdata
        if event == 'key_press_event': print(event)
        global coords, window, masks
        
        if len(window)%2 == 0: window = []
        window.append(ix)
        if len(window)%2 == 0:
            window = sorted(window)
            print(window)
            coords.append(window)
            m = plt.axvspan(window[0],window[1], alpha=0.3, color='grey')
            masks.append(m)
            plt.pause(0.05)
    if event.button == 3: 
        ix = event.xdata
        ToRemove = []
        for i,coord in enumerate(coords):
            if coord[0] <= ix and coord[1] >= ix:
                ToRemove.append([masks[i],coord])
        if len(ToRemove) > 0:
            for rem in ToRemove:
                rem[0].remove()
                masks.remove(rem[0])
                coords.remove(rem[1])
    return coords
def out(event):
    '''
    Function to leave the interavtive mode by pressing any key

    Returns
    -------
    None.

    '''
    fig.canvas.mpl_disconnect(cid)
    return None

def mask_func(wave,coords):
    '''
    Function to create the mask over the wavelength using the coordenates (output
    from "onclick")
    '''
    mask = False
    for coord in coords:
        mask += (wave >= coord[0]) & (wave <= coord[1])
    mask = np.invert(mask)
    return mask
    
##### Functions definition #####


# Test 
if __name__ == '__main__':
    ### Plotting observational spectrum ###
    star = sp.spec('HD49798')
    
    fig, ax = plt.subplots()
    
    ax.plot(star.waveflux()[0],star.waveflux()[1])
    ### Plotting observational spectrum ###
    
    #### Interactive selecion of the lines to be masked ####
    '''mask
    If pressing any key while in the figure, the interactive version is closed (not
    the figure)
    '''
    cid = fig.canvas.mpl_connect('button_press_event', onclick) # Select windows
    fig.canvas.mpl_connect('key_press_event', out) # Close interactive
    #### Interactive selecion of the lines to be masked ####
    print('\n')
#%% Run after finishing the interactive mode
#if __name__ == '__main__':
    mask_n = mask_func(star.waveflux()[0],coords)
    
    ax.plot(star.waveflux()[0][mask_n],star.waveflux()[1][mask_n])
