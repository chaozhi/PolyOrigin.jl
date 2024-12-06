function brentMin(f::Function,lowbound,upbound,xstart,precisiongoal,accuracygoal,maxiter)
    cgold = 0.381966 #(3-Sqrt(5))/2
    a = min(lowbound, upbound)
    b = max(lowbound, upbound)
    if isnothing(xstart) || (!(a<=xstart<=b))
        x = a + cgold*(b - a)
    else
        x = xstart
    end
    v = w = a + cgold*(b - a)
    e = 0
    d = 0 # required for compilation
    fv = fw = fx = f(x)
    his=zeros(Float64,maxiter+1,3)
    his[1,:] = [0,x,fx]
    for iter = 1:maxiter
        xm = 0.5*(a + b)
        tol1 = 10.0^(-precisiongoal)*abs(x) + 10.0^(-accuracygoal)
        tol2 = 2*tol1
        # abs(x - xm) <= tol2 - 0.5*(b - a) <==> max(x-a,b-x)<=tol2
        if abs(x - xm) <= tol2 - 0.5*(b - a)
            his = his[1:iter,:]
            break
        end
        p = q = r = 0.
        if abs(e) > tol1
            #Parabolic fit from x, v, and w
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2*(q - r)
            #result in proposal step size -p/q with always q\[GreaterEqual]0
            if q > 0
                p = -p
            else
                q = -q
            end
            # r denotes the movement of the step before last
            r = e
            #e denotes the movement of the last step
            e = d
        end
        if abs(p) < abs(0.5*q*r) && p > q*(a - x) && p < q*(b - x)
            #Parabolic interpolation step
            d = p/q
            u = x + d
            # f must not be evalued too close to bracketing ends a or b
            if u - a < tol2 || b - u < tol2
                d = x<xm ? tol1 : -tol1
            end
        else
            # Golden section step
            e = (x<xm ? b : a)-x
            d = cgold*e
        end
        # f must not be evalued too close to x
        u = x+ (abs(d) >= tol1 ? d : (d>0 ? tol1 : -tol1))
        fu = f(u)
        # housekeeping
        if fu <= fx
            if u < x
                b = x
            else
                a = x
            end
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
        else
            if u < x
                a = u
            else
                b = u
            end
            if fu <= fw || w == x
                v = w
                fv = fw
                w = u
                fw = fu
            else
                if fu <= fv || v == x || v == w
                    v = u
                    fv = fu
                end
            end
        end
        his[iter+1,:] = [iter,x,fx]
    end
    his = vcat(["iteration" "x" "f(x)"],his)
    x,fx,his
end

function brentMin(f::Function,lowbound,upbound,xstart=nothing;
    precisiongoal::Real=6,accuracygoal::Real=6,maxiter::Integer=100)
    brentMin(f,lowbound,upbound,xstart,precisiongoal,accuracygoal,maxiter)
end

function brentMax(f::Function,lowbound,upbound,xstart=nothing;
    precisiongoal::Real=6,accuracygoal::Real=6,maxiter::Integer=100)
    f2(x) =-f(x)
    (x,fx,his) = brentMin(f2,lowbound,upbound,xstart,precisiongoal,accuracygoal,maxiter)
    fx*=-1
    his[2:end,3]*=-1
    x,fx,his
end
