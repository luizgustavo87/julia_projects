# UASB 2D
using Gridap
using SparseArrays
using DifferentialEquations
using BenchmarkTools

# Parâmetros do problema
Re = 300.0  # Número de Reynolds
Sc = 1.0    # Número de Schmidt
Fr = 0.01    # Número de Froude
g = [0.0, -1.0]  # Vetor gravidade
Kdis=0.5

# Define domain and mesh
n = 30
Ly = 3
domain = (0, 1, 0, Ly)
partition = (n, 3*n)
model = CartesianDiscreteModel(domain, partition; isperiodic=(true, false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "diri1", [6])
add_tag_from_tags!(labels, "diri0", [1, 2, 3, 4, 5, 7, 8])
add_tag_from_tags!(labels, "diri00", [1, 2, 3, 4, 7, 8])

add_tag_from_tags!(labels, "top", [6])
add_tag_from_tags!(labels, "botton", [5])
add_tag_from_tags!(labels, "right", [8])
add_tag_from_tags!(labels, "left", [7])
add_tag_from_tags!(labels, "corner", [1, 3])

D = 2
order = 2 #ordem do polinômio que aprox a vel.
reffeᵤ = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, labels=labels, dirichlet_tags=["botton"]) #campo de velocidade

reffeₚ = ReferenceFE(lagrangian, Float64, order-1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, dirichlet_tags=["top"]) #campo de pressão e function de pressões

reffec1 = ReferenceFE(lagrangian, Float64, order)
R1 = TestFESpace(model, reffec1, conformity=:H1, labels=labels, dirichlet_tags=["botton"]) #campo de concentração 1

reffec2 = ReferenceFE(lagrangian, Float64, order)
R2 = TestFESpace(model, reffec2, conformity=:H1, labels=labels, dirichlet_tags=["botton"]) #campo de concentração 2

uD0 = VectorValue(0, 0) 
uD1 = VectorValue(0, 1)
cD1 = 1.0

#################################################
#  Funções de Condições Iniciais e Contorno (Input the initial conditions and boundary conditions)

g0(x, t::Real) = uD0 #Define uma condição de contorno de velocidade constante uD0 para todos os pontos x no tempo t. VectorValue(0, 0)
g0(t::Real) = x -> g0(x, t) #Cria uma função que retorna a condição de contorno de velocidade para qualquer ponto x no tempo t.

g1(x, t::Real) = VectorValue(0, 10.0 * (1.0 - sqrt(x[2] / Ly)) * (((0.5 * (1.0 - (cos(3 * 2 * pi * x[1]))) + 0.0 * rand()) > 0.99 ? 1.0 : 0.0)) + sqrt(x[2] / Ly)) #O 3 é numero de jatos
g1(t::Real) = x -> g1(x, t)
p1(t::Real) = x -> 0.0 #0.0 na direção x
#c0_1(x, t::Real) = 1.0
#c0_1(t::Real) = x -> c0_1(x, t)
c1(x, t::Real) = 1.0 #colocar o valor real ou 1 mesmo.
c1(t::Real) = x -> c1(x, t)
#c0_2(x, t::Real) = 0.0
#c0_2(t::Real) = x -> c0(x, t)
c2(x, t::Real) = 0.0
c2(t::Real) = x -> c2(x, t)

#c0(x, t::Real) = 0.0
#c0(t::Real) = x -> c0(x, t)

#Transient Trial FE Space
U = TransientTrialFESpace(V, [g1]) #cond de contorno de vel
P = TransientTrialFESpace(Q)
C1 = TransientTrialFESpace(R1, [c1])
C2 = TransientTrialFESpace(R2, [c2])

Y = MultiFieldFESpace([V, Q, R1, R2])
X = TransientMultiFieldFESpace([U, P, C1, C2])

degree = 2*order #Confirmar com o Prof. Norberto se é =order ou 2*order (CONFIRMADO!!)
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ, degree) #usar p/ todos 

Γₕ = BoundaryTriangulation(model, tags="botton") #todo o botton
dΓₕ = Measure(Γₕ, degree)
n_Γₕ = get_normal_vector(Γₕ)

s1 = VectorValue(0.0, 1.1) #veloc de sendiment (queda) - maior que a do escoamento e tende a cair 1. Onde o jato entra ela sobe.
s2 = VectorValue(0.0, 1.1)
c_t = 1.0   
beta = 0.1 
conv(u, ∇u) = (∇u') ⋅ u
buss(v,c1,c2) = beta * v[2] * (c1+c2) #verificar se é const ou var
wwc1(u) = u - s1  #wwc é o up_i no mini artigo! 
wwc2(u) = u - s2
F1(c1)=-Kdis*c1
F2(c1)=Kdis*c1
#convc1(u, ∇c1) = (∇c1) ⋅ u
#convc2(u, ∇c2) = (∇c2) ⋅ u

a((u, p, c1, c2), (v, q, r1, r2)) = ∫( 1/Re*(∇(v)⊙∇(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ + 
                          ∫( 1/(Re*Sc)*(∇(r1)⊙∇(c1)) )dΩ + 
                          ∫( 1/(Re*Sc)*(∇(r2)⊙∇(c2)) )dΩ +
                          ∫( v ⋅ (conv∘(u, ∇(u))) )dΩ + 
                          ∫( buss∘(v,c1,c2) )dΩ +
                          ∫( r1 ⋅ inner(wwc1(u), ∇(c1)) )dΩ + 
                          ∫( r2 ⋅ inner(wwc2(u), ∇(c2)) )dΩ - 
                          ∫((F1∘(c1))*r1)dΩ - 
                          ∫((F2∘(c1))*r2)dΩ +
                          ∫( (v ⋅ n_Γₕ) * (beta * c1) )dΓₕ +
                          ∫( (v ⋅ n_Γₕ) * (beta * c2) )dΓₕ

at(t, (u, p, c1, c2), (v, q, r1, r2)) = ∫( ∂t(u)⋅v )dΩ +
                                            ∫( ∂t(c1)⋅r1 )dΩ +
                                            ∫( ∂t(c2)⋅r2 )dΩ

#cc(u, (c1, c2, r1, r2)) = ∫( r1 ⊙ (convc1∘(u, ∇(c1))) )dΩ +
                       #∫( r2 ⊙ (convc2∘(u, ∇(c2))) )dΩ

#dc(u, du, v) = ∫( v ⊙ (dconv∘(du, ∇(du), u, ∇(u))) )dΩ

res(t, (u, p, c1, c2), (v, q, r1, r2)) = at(t, (u, p, c1, c2), (v, q, r1, r2)) + a((u, p, c1, c2), (v, q, r1, r2))
#jac((u,p, c1, c2),(du,dp,dc1,dc2),(v,q,r1,r2)) = a((du,dp,dc1,dc2),(v,q)) + dc(u,du,v) 
                                          
op = TransientFEOperator(res, X, Y)

using LineSearches: BackTracking
nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=30)

CFL = 1.0/10
Δt = CFL*1/n
θ = 0.5
t₀ = 0.0
n_iteracoes = 30
T = n_iteracoes*Δt 
ode_solver = ThetaMethod(nls,Δt,θ)



@info "Calculando campos de velocidade e pressão..."


# ThetaMethod is not directly compatible with DifferentialEquations.jl
# We will redefine the ODE function and solve the problem using DifferentialEquations.jl
# function ODEfun(du, u, p, t)
#     du .= zero(u)
#     residual = op(t, u)
#     for i in 1:length(du)
#         du[i] = residual[i]
#     end
# end

# Define initial condition U₀ as a vector
U₀ = interpolate_everywhere([g1(0), p1(0), c2(0)], X(0.0))
uhs = Gridap.solve(ode_solver, op, U₀, t₀, T) 
# Define the ODE problem
# prob = ODEProblem(ODEfun, U₀, (t₀, T))

# Solve the problem
# sol = OrdinaryDiffEq.solve(prob, Tsit5())
# sol = solve(prob, Tsit5())

# Export results
output_dir = "D:\\Mestrado\\2024-1\\Seminários\\Mariana\\desintegration_results"
# for (i, uₕ) in enumerate(sol.u)
#     uh, ph, ch = uₕ
#     if mod(i, 10) == 0
#         writevtk(Ωₕ, "$output_dir/results2_$(div(i, 10)).vtu", cellfields=["uh" => uh, "ph" => ph, "ch" => ch])
#     end
# end

it=1
  for (uhit,t) in uhs
    global it
	
    @info "Iteração atual: $it"
    uh, ph, ch1, ch2 = uhit
    writevtk(Ωₕ, "$output_dir/desintegration_$it.vtu", cellfields=["uh" => uh, "ph" => ph, "ch1" => ch1, "ch2" => ch2])

    it += 1
  end