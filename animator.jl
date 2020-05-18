using Plots
using CSV

DATAFILE = "runs/ManouvreTest.csv"
SPEED = 100

data = CSV.File(DATAFILE)

plt = plot3d(
    1,
    xlim = (-30000, 30000),
    ylim = (-30000, 30000),
    zlim = (-30000, 30000),
    title = "Anim",
    grid = nothing,
    axes = false,
    ticks = false,
    )
    
anim = @animate for row in data[1:SPEED:end]
    push!(plt, row.x, row.y, row.z)
end every 10
# surface(x,y,z, color="red")

gif(anim, "animations/ManouvreTest.gif")