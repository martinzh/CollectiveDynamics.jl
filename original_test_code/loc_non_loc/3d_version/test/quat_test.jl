using Plots, Quaternions
gr()
pyplot()

q0 = Quaternion([1., 0., 0.])

ax = [-1., 1., 1.]

x = Float64[]
y = Float64[]
z = Float64[]

x_1 = Float64[]
y_1 = Float64[]
z_1 = Float64[]

push!(x, q0.v1)
push!(y, q0.v2)
push!(z, q0.v3)

push!(x_1, q0.v1)
push!(y_1, q0.v2)
push!(z_1, q0.v3)

for i in 0:0.05:π

    qR = qrotation(cross([q0.v1, q0.v2, q0.v3], ax), i) * q0
    qR_1 = q0 * qrotation(cross([q0.v1, q0.v2, q0.v3], ax), i)

    push!(x, qR.v1)
    push!(y, qR.v2)
    push!(z, qR.v3)

    push!(x_1, qR_1.v1)
    push!(y_1, qR_1.v2)
    push!(z_1, qR_1.v3)

end

scatter(x, y, z)

scatter(hcat(x, x_1), hcat(y, y_1), hcat(z, z_1))


gui()
