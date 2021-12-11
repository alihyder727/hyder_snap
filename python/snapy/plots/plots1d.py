from numpy import mean, std, log

def PlotProfile1(ax, var, yaxis, **args):
  if len(var.shape) > 1:
    a = list(range(len(var.shape)))
    a.remove(1)
    var_avg = mean(var, axis = tuple(a))
    var_std = std(var, axis = tuple(a))
    h = ax.plot(var_avg, yaxis, **args)
    ax.fill_betweenx(yaxis, var_avg - var_std, var_avg + var_std,
      facecolor = h[0].get_color(), edgecolor = None, alpha = 0.2)
  else:
    ax.plot(var, yaxis, **args)

