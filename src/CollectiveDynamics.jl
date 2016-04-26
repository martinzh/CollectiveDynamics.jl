# Modelo movimiento colectivo
# Paralelo con SharedArrays

module CollectiveDynamics

export Param

#///////////////////////////#

type Param
    m::Int64
    r0::Float64
    eta::Float64
    w::Float64
    bound::Float64

    # Inicializa todo a ceros
    Param() = new(0, 0.0, 0.0, 0.0, 0.0)
end

#///////////////////////////#

function init_pos(N::Int64, double w, double* pos, gsl_rng* r) {
    int i;
    double u;
    for (i = 0; i < 2*N; i++) {
        u = randNum(w, r);
        pos[i] = u;
    }
}

#///////////////////////////#

function init_vels(N::Int64, double* v, double v0, gsl_rng* r) {
    int i;
    float vx, vy, normV;

    for (i = 0; i < N; i++) {
        vx = randNum(1.0, r);
        vy = randNum(1.0, r);

        normV = sqrt(vx*vx + vy*vy);

        v[2*i]   = v0 * (vx / normV);
        v[2*i+1] = v0 * (vy / normV);
    }
}


#///////////////////////////#

end
