using LinearAlgebra, Plots

cd("My PATH")
wdir = pwd();

# Import function that does the iteration/interpolation steps for us
include("model_functions.jl")


# Part I - Set up numerical calculation

# Set parameters
β     = 0.95;           # Discount factor (Mark Whatney's impatience)
α     = 1/3;            # Output elasticity of capital (raw potatoes)
δ     = 0.5;              # Depreciation: Share of capital (raw potatoes) that deteriorates
                        # every period. Set between 1 (full) or greater 0.

N     = 1000;           # Choose number of grid points
Kmin  = 0.01;            # lowest point on capital grid
Kmax  = 1;            # highest point on capital grid
CapC  = exp(-1e16);     # force lower bound on consumption policy

Kgrid = range(Kmin, stop = Kmax, length = N);

## Part II - Calculate VF numerically and find minimum using Brent method

# Do value-function iteration for given K on Kgrid
d(K) = VF_solver(K, Kgrid, α, δ, β, CapC)

Kstar, iter = Brent(d, Kmin, Kmax)

# Find aggregate number of produced potatoes, consumable potatoes and necessary
# investment (new planting of potatoes) for next period.

Ystar  = Kstar.^(α);
Cstar  = Ystar-δ*Kstar;
I_star = δ*Kstar;

## Part III: Compare analytical and numerical solution (when δ = 1, i.e. full depreciation).

Kstar_an  = (α/(1/β-(1-δ))).^(1/(1-α));
Ystar_an  = Kstar_an.^(α);
Cstar_an  = Ystar_an-δ*Kstar_an;
I_star_an = δ*Kstar_an;

δ

a             = 1/(1-β)*(log(1-α*β)+(α*β*log(α*β)/(1-α*β)));
phi           = α/(1-β*α);

VF_analytical = phi*log.(Kgrid).+a; # Value function
Kpol_an       = β*phi/(1+β*phi).*Kgrid.^(α) # Capital policy
Cpol_an       = Kgrid.^(α).+(1-δ).*Kgrid.-Kpol_an # Consumption Policy

# Evaluate Results again to plot numerical solution for Value- and
# Policy-Function against analytical result (when when δ = 1, i.e. full depreciation).

excess_supply, Kstar, Vnew, pol, Kpol, Cpol = VF_solver(Kstar, Kgrid, α, δ, β, CapC)

## Part IV: Plot solution

pyplot()

# Value Function

# analytical part
plot(Kgrid,VF_analytical,linewidth = 2, xlabel = "Asset Grid", ylabel = "Value Function", label = "VF Analytical")
# add numerical part
plt1 = plot!(Kgrid,Vnew,linewidth = 2,  xlabel = "Asset Grid", ylabel = "Value Function", label = "Numerical", title = "Comparing VFIs")

# Policy Function

plot(Kgrid,Cpol_an,     linewidth = 2, xlabel = "Asset Grid", ylabel = "Consumption Policy", label = "Policy Analytical")
plt2 = plot!(Kgrid,Cpol,linewidth = 2, xlabel = "Asset Grid", ylabel = "Consumption Policy", label = "Numerical", title = "Comparing Policies")

# Combine in one Graph

plt3 = plot(plt1,plt2, layout = (1,2),legend = true, dpi = 300)

savepath = wdir*"/plots/combined_"*string(N)*".png"
savefig(savepath)

(δ/α).^(1/(α-1))

c_steady = Kgrid.^(α).-δ.*Kgrid

plt4 = plot(Kgrid,c_steady)
plot!([Kstar_an], seriestype="vline", label="", legend = false, ylabel = "Potato Consumption", xlabel = "Raw Potato Grid", dpi = 300, title = "Mark's Optimal Consumption Policy", grid =true, gridlinewidth=1, tick_direction =:out, foreground_color_legend = nothing)
annotate!(Kstar_an, 0, text("\$ K^\\ast \$", :bottom));
display(plt4)

savepath = wdir*"/plots/optimal_consumption_"*string(N)*".png"
savefig(savepath)


# Reproduce the declining potato consumption every period:

days = 4*365+1
β = 0.005;
n_potatoes_per_day = 10;

cons_over_time     = zeros(days,1)
cons_over_time[1]  = 1000;

for i_cons = 2:days
    cons_over_time[i_cons] = (1-β).*cons_over_time[i_cons-1]
end

plot(1:days,cons_over_time, xlabel = "No. of Days until Rescue", ylabel = "Number of Potatoes", legend = false, title = "Mark's time until he runs out of Potatoes",grid =true, gridlinewidth=1, tick_direction =:out, foreground_color_legend = nothing, tickfont=font(10), ytickfont=font(10), guidefont=font(10), legendfont=font(10), formatter = :plain, dpi = 300, size = (800, 300))

savepath = wdir*"/plots/runout_potatoes"*".png"
savefig(savepath)
