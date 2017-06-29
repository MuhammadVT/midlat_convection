def genCmap(scale, num_seg=8, colors='lasse'):
    """Generates a colormap and returns the necessary components to use it
    Parameters
    ----------
    scale : list
        a list with the [min,max] values of the color scale
    colors : Optional[str]
        a string indicating which colorbar to use, valid inputs are 
        'lasse', 'mr'.  default = 'lasse'
    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        the colormap generated.  This then gets passed to the mpl plotting
        function (e.g. scatter, plot, LineCollection, etc.)
    norm : matplotlib.colors.BoundaryNorm
        the colormap index.  This then gets passed to the mpl plotting
        function (e.g. scatter, plot, LineCollection, etc.)
    bounds : list
        the boundaries of each of the colormap segments.  This can be used
        to manually label the colorbar, for example.
    """

    import matplotlib,numpy
    import matplotlib.pyplot as plot

    #the MPL colormaps we will be using

    cmj = matplotlib.cm.jet
    cmpr = matplotlib.cm.prism

    #check for a velocity plot

    #check for what color scale we want to use
    if(colors == 'mr'):
        #define our discrete colorbar
        if num_seg == 6:
        cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                     cmpr(.11), cmpr(.1),
                                                     cmpr(.175), cmpr(.158),
                                                     cmj(.32), cmj(.37)])

        elif num_seg == 7:
        cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                     cmpr(.11), cmpr(.1),
                                                     cmpr(.175), cmpr(.158),
                                                     cmj(.32), cmj(.37)])


        elif num_seg == 8:
        cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                     cmpr(.11), cmpr(.1),
                                                     cmpr(.175), cmpr(.158),
                                                     cmj(.32), cmj(.37)])

    else:
        #define our discrete colorbar
        if num_seg == 6:
        cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8),
                                                 cmj(.7), cmj(.65),
                                                 cmpr(.142), cmj(.45),
                                                     cmj(.3), cmj(.1)])

        elif num_seg == 7:
        cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8),
                                                 cmj(.7), cmj(.65),
                                                 cmpr(.142), cmj(.45),
                                                     cmj(.3), cmj(.1)])
        elif num_seg == 8:
        cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8),
                                                 cmj(.7), cmj(.65),
                                                 cmpr(.142), cmj(.45),
                                                     cmj(.3), cmj(.1)])


    #define the boundaries for color assignments
    bounds = numpy.round(numpy.linspace(scale[0],scale[1],num_seg))
    bounds = numpy.append(bounds,10000.)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)


    # If its a non-velocity plot, check what color scale we want to use
    if(colors == 'mr'):
        #define our discrete colorbar
        cmap = matplotlib.colors.ListedColormap([cmpr(.175), cmpr(.158),
                                                 cmj(.32), cmj(.37),
                                                 cmpr(.142), cmpr(.13),
                                                 cmpr(.11), cmpr(.10)])
    else:
        #define our discrete colorbar
        cmap = matplotlib.colors.ListedColormap([cmj(.1), cmj(.3), cmj(.45),
                                                 cmpr(.142), cmj(.65),
                                                 cmj(.7), cmj(.8), cmj(.9)])

    #define the boundaries for color assignments
    bounds = numpy.round(numpy.linspace(scale[0],scale[1],8))
    bounds = numpy.append(bounds,10000.)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    cmap.set_bad('w',1.0)
    cmap.set_over('w',1.0)
    cmap.set_under('.6',1.0)

    return cmap,norm,bounds

def addColorbar(mappable, ax):
    """ Append colorbar to axes

    Parameters
    ----------
    mappable :
        a mappable object
    ax :
        an axes object

    Returns
    -------
    cbax :
        colorbar axes object

    Notes
    -----
    This is mostly useful for axes created with :func:`curvedEarthAxes`.

    written by Sebastien, 2013-04

    """
    from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size
    import matplotlib.pyplot as plt 

    fig1 = ax.get_figure()
    divider = SubplotDivider(fig1, *ax.get_geometry(), aspect=True)

    # axes for colorbar
    cbax = LocatableAxes(fig1, divider.get_position())

    h = [Size.AxesX(ax), # main axes
         Size.Fixed(0.1), # padding
         Size.Fixed(0.2)] # colorbar
    v = [Size.AxesY(ax)]

    _ = divider.set_horizontal(h)
    _ = divider.set_vertical(v)

    _ = ax.set_axes_locator(divider.new_locator(nx=0, ny=0))
    _ = cbax.set_axes_locator(divider.new_locator(nx=2, ny=0))

    _ = fig1.add_axes(cbax)

    _ = cbax.axis["left"].toggle(all=False)
    _ = cbax.axis["top"].toggle(all=False)
    _ = cbax.axis["bottom"].toggle(all=False)
    _ = cbax.axis["right"].toggle(ticklabels=True, label=True)

    _ = plt.colorbar(mappable, cax=cbax)

    return cbax

