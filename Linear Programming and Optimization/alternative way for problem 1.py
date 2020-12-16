def sde(x):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    return x.std()

def parset(x, method = 'Mean'):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    if method.upper() == 'MEAN':
        return x.mean()
    if method.upper() == 'MEDIAN':
        return x.median()

def bsample(x, size = None):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    if size == None:
        size = x.shape[0]
    return x.sample(n = size, replace = True)


def bootpar(x, method = 'Mean', size = None):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    y = bsample(x, size)
    m = parset(y, method)
    return m
    
def repeat_fun(f, times, x, method = 'Mean', size = None):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return    
    return pd.DataFrame([f(x, method = method, size = size) for _ in range(times)])

def boot_bias_sde(x, boottimes = 10, method = 'Mean', size = None):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    a = repeat_fun(bootpar, times = boottimes, x = x, method  = method, size = size)
    if method.upper() == 'MEAN':
        b = x.mean()
    if method.upper() == 'MEDIAN':
        b = x.median()
    botse = a.std()
    bias = a.mean() - b
    return {"Original": b[0], "Bias": bias[0], "Standard Error": botse[0]}