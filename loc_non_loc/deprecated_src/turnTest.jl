include("obj_lib.jl")
using PyPlot

function movePart(part::Bird,dt::Float64)
  part.pos += scale(part.vel,dt)
end

function resetPart(part::Bird)
  part.pos = [0.0,10.0]
  part.vel = [1.0,0.0]
end

function RotVec!(vec::Array{Float64,1},alpha::Float64)
  X = vec[1]*cos(alpha) - vec[2]*sin(alpha)
  Y = vec[1]*sin(alpha) + vec[2]*cos(alpha)

  vec[1] = X
  vec[2] = Y
end

function makeTurn(vel::Array{Float64},rate::Float64)
  ang = (0.5*pi)/rate
  RotVec!(vel,ang)
end

function makeTray(part::Bird,t::Int64)
  tray = zeros(t,2)
  resetPart(part)

  for i in 1:t
    movePart(part,1.0)
    tray[i,1] = part.pos[1]
    tray[i,2] = part.pos[2]
  end
  return tray
end

function makePertTray(part::Bird,t::Int64,rate::Float64
                      ,τ::Float64,velPert::Array{Float64,1})
#function makePertTray(part::Bird,t::Int64,rate::Float64,τ::Float64)

  resetPart(part)
  println("reset:",part.vel)

  tray = zeros(t,2)
  velPert = zeros(2)

  for i in 1:t

    if i == τ
      println("antes de copiar")
      println(velPert,"\t",part.vel)
      copy!(velPert,part.vel)
#      velPert = copy(part.vel)
      println("despues de copiar")
      println(velPert,"\t",part.vel)
    end

    if i>= τ && i < τ+rate
      makeTurn(velPert,rate)
      #part.vel = copy(velPert)
      copy!(part.vel,velPert)
    end

    movePart(part,1.0)

    tray[i,1] = part.pos[1]
    tray[i,2] = part.pos[2]
  end
  return tray
end

part = Bird(zeros(2),[1.0,0.0],[0])
velPert = Array(Float64,2)

tray = makeTray(part,50)
trayPert = makePertTray(part,50,35.0,25.0,velPert)

plt.plot(tray[:,1],tray[:,2],".-")
plt.plot(trayPert[:,1],trayPert[:,2],".-")
plt.clf()

tray
trayPert

