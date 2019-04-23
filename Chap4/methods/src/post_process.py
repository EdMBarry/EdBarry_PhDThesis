
def load_result(filename):
    ''' Loads in the data from the 2D simulations.
    filename: Input the relative path of the csv file you want to import.
    
    We drop the uninteresting parts of the results.
    '''
    import pandas as pd
    
    result = pd.read_csv(filename).drop(["ALambda_0", "D_X","SAC",
                                  "SI","SIC","SIP","SLambda_0",
                                  "SS","U:0","U:1","U:2","Points:0","Points:1",
                                  "Points:2","vtkValidPointMask","p","p_rgh",
                                  "XI_VS","XPB_VS","XS_VS"],
                                   axis=1)
    return result


def plot_subplot(dataset0, dataset1, dataset2, dataset3):
    ''' Plot the data for each timestep of interest
    datasetN: define 4 datasets for the subplot
    '''
    files = [dataset0, dataset1, dataset2, dataset3]
    days = ["Time: 1.0 d", "Time: 4.0 d", "Time: 7.0 d", "Time: 10.0 d"]
    
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import matplotlib as mpl
    
    # Set the font parameters and linewidths
    mpl.rcParams['lines.linewidth'] = 3
  #  mpl.rcParams.update({'figure.autolayout': True})
    set_color="#444444"
    plt.style.use("seaborn-whitegrid")
    font = {'family':'DejaVu Sans','weight':'bold','size':16}
    plt.rc('font',**font)
    
    # Start plotting the figures
    fig = plt.figure(figsize=(10,7.5), frameon=False, dpi=200)

    for i in range(4):
        ax = fig.add_subplot(2,2,i+1)
        ax.plot(files[i]['arc_length']*1000, files[i]['XPB_alpha'],
                 color="#ff0000",
                 label=r"$\mathbf{\alpha_{XPB}}$")
        
        ax.plot(files[i]['arc_length']*1000, files[i]['XS_alpha'],
                 color="#1aff1a",
                 label=r"$\mathbf{\alpha_{XS}}$")
        
        ax.plot(files[i]['arc_length']*1000, files[i]['XI_alpha'],
                 color="#996600",
                 label=r"$\mathbf{\alpha_{XI}}$")
        
        ax.plot(files[i]['arc_length']*1000, files[i]['alpha.water'],
                 color="#1a75ff",
                 label=r"$\mathbf{\alpha_{water}}$")
        
        plt.title(days[i], fontsize=12, fontweight="bold")
        plt.xlim([0.0,3.0])
        plt.ylim([-0.01,1.01])

        plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
        plt.locator_params(axis='y', nbins=5)    
        plt.gca().xaxis.get_majorticklabels()[0].set_x(-0.1)


    plt.text(-2, -0.2, "Distance from substratum (mm)")
    plt.text(-4.5, 1.2, "Volume fraction (-)", va='center', rotation='vertical')
    handles, labels = ax.get_legend_handles_labels()

    fig.tight_layout()
    plt.legend(handles, labels, fontsize=12, loc=(-0.7, -0.3), ncol=4)