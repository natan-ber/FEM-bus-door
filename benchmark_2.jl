using Gridap, GridapGmsh, LinearAlgebra, SparseArrays, Preconditioners, Gridap.Geometry, IterativeSolvers, Random, Test, BenchmarkTools
model = GmshDiscreteModel("C:/.../bus_door.msh")
# writevtk(model,"model")

labels = get_face_labeling(model)
dimension = 3
tags = get_face_tag(labels,dimension)

const glass_tag = get_tag_from_name(labels,"glass")

order = 1

reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V0 = TestFESpace(model,reffe;
  conformity=:H1,
  dirichlet_tags=["boundary1","boundary2","boundary3"],
  dirichlet_masks=[(true,true,true), (true,true,true), (false,false,true)])
  
g1(x) = VectorValue(0.0,0.0,0.0)
g2(x) = VectorValue(0.0,0.0,0.0)
g3(x) = VectorValue(0.0,0.0,-0.005)

U = TrialFESpace(V0,[g1,g2,g3])

function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

const E_glass = 100e9
const ν_glass = 0.37
const (λ_glass,μ_glass) = lame_parameters(E_glass,ν_glass)

const E_alu = 69.0e9
const ν_alu = 0.33
const (λ_alu,μ_alu) = lame_parameters(E_alu,ν_alu)

function σ_bimat(ε,tag)
  if tag == glass_tag
    return λ_glass*tr(ε)*one(ε) + 2*μ_glass*ε
  else
    return λ_alu*tr(ε)*one(ε) + 2*μ_alu*ε
  end
end

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags)) )*dΩ
l(v) = 0

op = AffineFEOperator(a,l,U,V0);
A = get_matrix(op);
b = get_vector(op);
p = AMGPreconditioner{SmoothedAggregation}(A);
@benchmark x = cg(A,b,verbose=true,Pl=p)
x = cg(A,b,verbose=true,Pl=p);
uh = FEFunction(U,x);

writevtk(Ω,"results",cellfields=
 ["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ_bimat∘(ε(uh),tags)])
