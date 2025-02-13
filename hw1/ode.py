def rk4(f,y0,t0,h,N):
    y = [y0]
    t = [t0]
    for i in range(N):
        k1 = f(t[-1],y[-1])
        k2 = f(t[-1]+h/2,y[-1]+h/2*k1)
        k3 = f(t[-1]+h/2,y[-1]+h/2*k2)
        k4 = f(t[-1]+h,y[-1]+h*k3)
        y.append(y[-1] + (h/6)*(k1 + 2*(k2 + k3) + k4))
        t.append(t[-1] + h)
    return t,y