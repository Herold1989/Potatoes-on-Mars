function VFI(V, u, β, dist)

    while dist > 1e-10

    Vlong      = repeat(V,1,length(V))';
    Vnew,pol   = findmax(real(u.+β.*Vlong),dims = 2);

    dist = maximum(broadcast(abs, (Vnew-V))); # calculate distance
    V    = Vnew; # update VF

    end

    Vlong      = repeat(V,1,length(V))';
    Vnew,pol   = findmax(real(u.+β.*Vlong),dims = 2);

    return Vnew, pol
end

function mylinearinterpolate(xgrd,ygrd, xeval)

    ind = zeros(Int,length(xeval))
    n_xgrd = length(xgrd)
    yeval= zeros(eltype(ygrd),length(xeval))
    #ind = zeros(Int,length(xeval))
    n_xgrd = length(xgrd)
    @views for i in eachindex(xeval)
        xi = xeval[i]
        if xi .> xgrd[n_xgrd-1]
            iL = n_xgrd - 1
        elseif xi .< xgrd[2]
            iL = 1
        else
            iL = locate(xi, xgrd)
        end
        iR = iL+1
        xL = xgrd[iL]
        wR = (xi .- xL)./ (xgrd[iR] .- xL)
        wL = 1.0-wR
        yeval[i] = wL .* ygrd[iL] .+ wR .* ygrd[iR]
    end

    return yeval
end


#-------------------------------------------------------------------------------
## BRENT'S METHOD ##
#-------------------------------------------------------------------------------
function Brent(f, a, b)

    tol = 1e-14;
    # Implementation of Brent's method to find a root of a function (as on wikipedia)
    fa = f(a)[1]; fb = f(b)[1]
    if fa * fb > 0
        error("f[a] and f[b] should not have different signs!")
    end

    c = a;
    fc=fa;   # at the beginning: c = a
    c = a;
    d = b - a;
    e = d

    iter = 0
    maxiter = 1000

    while iter < maxiter
        iter += 1

        if fb * fc > 0
            c = a; fc = fa; d = b - a; e = d
        end

        if abs(fc) < abs(fb)
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
        end
        tol = 2.0 * eps() * abs(b) + tol; m = (c - b) / 2.0; #Toleranz

        if abs(m)>tol && abs(fb)>0 #Verfahren muss noch durchgeführt werden
            if abs(e)<tol || abs(fa)<=abs(fb)
                d=m; e=m
            else
                s=fb/fa
                if a==c
                    p=2*m*s; q=1-s
                else
                    q=fa/fc; r=fb/fc
                    p=s*(2*m*q*(q-r)-(b-a)*(r-1))
                    q=(q-1)*(r-1)*(s-1)
                end
                if p>0
                    q=-q
                else
                    p=-p
                end
                s=e; e=d
                if  2*p<3*m*q-abs(tol*q)  && (p<abs(s*q/2))
                    d=p/q
                else
                    d=m; e=m
                end
            end
            a=b; fa=fb
            if abs(d)>tol
                b=b+d
            else
                if m>0
                    b=b+tol
                else
                    b=b-tol
                end
            end
        else
            break
        end
        fb = f(b)[1]
    end
    return b, iter
end

function VF_solver(K_guess,Kgrid,α,δ,β,CapC)

        V = log.(Kgrid) # first (primitive) guess of VF
        N = length(Kgrid)

    # 2D - Consumption grid: Different consumption choices for different capital states
    C_aux = repeat((Kgrid.^α+(1-δ)*Kgrid),1,N)-repeat(Kgrid',N,1);
    # Define a utility function
                 u = log.(Complex.(C_aux));
    # Force lower bound on utility from zero consumption
    u[C_aux.<0.0] .= log(CapC);
                 u = real(u); # get rid of imaginary part
              dist = 5;
         Vnew, pol = VFI(V, u, β, dist)
                 x = K_guess;
          GridData = Kgrid[getindex.(pol,2)];
           Savings = mylinearinterpolate(Kgrid, GridData, x)
     excess_supply = K_guess - Savings[1];
             Kpol  = Kgrid[getindex.(pol,2)]
             Cpol  = Kgrid.^(α).+(1-δ).*Kgrid.-Kpol

    return excess_supply, Savings, Vnew, pol, Kpol, Cpol

end

#-------------------------------------------------------------------------------
## Locate function ##
#-------------------------------------------------------------------------------
locate(x::Number, xx::AbstractVector) = exp_search(x, xx)
function bin_search(x::Number, xx::AbstractVector)
    # real(8), intent(in), dimension(N)		:: xx  ! lookup table
    # integer, intent(in)						:: N   ! no elements of lookup table
    # real(8), intent(in)						:: x   ! Value whose nearest neighbors are to be found
    # integer, intent(out)					:: j   ! returns j if xx(j)<x<xx(j+1),
    #												   ! 0 if value is to the left of grid,
    #												   ! N if value is to the right of grid

    N = length(xx)
    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        jl = 1
        ju = N
        while (ju-jl) != 1
            jm = div((ju+jl),2)

@inbounds   if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end

    return j
end

function exp_search(x::Number, xx::AbstractVector)
    # real(8), intent(in), dimension(N)		:: xx  ! lookup table
    # integer, intent(in)						:: N   ! no elements of lookup table
    # real(8), intent(in)						:: x   ! Value whose nearest neighbors are to be found
    # integer, intent(out)					:: j   ! returns j if xx(j)<x<xx(j+1),
    #												   ! 0 if value is to the left of grid,
    #												   ! N if value is to the right of grid
@inbounds begin
    N = length(xx)
    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        bound = 1
        while bound < N && x>xx[bound]
            bound *=2
        end
        jl = div(bound,2)
        ju = min(N,bound)
        while (ju-jl) != 1
            jm = div((ju+jl),2)

            if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end
end
    return j
end
