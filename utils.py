#Some useful functions

def Colornum_Mos_Metro(num):
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
    9:"gray",
    0:"lawngreen"}.get(num%10)  


def time2ind(t):
    #
    return np.where(neurotime >= t)[0][0]