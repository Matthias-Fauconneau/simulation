#![allow(incomplete_features,dead_code,non_snake_case)]
#![feature(const_generics,non_ascii_idents,fn_traits,unboxed_closures,trait_alias,const_compare_raw_pointers,box_syntax,option_unwrap_none)]
use std::cmp::max;
use framework::{core::{Zero,mask,cb,sqrt}, vector::{xy,uint2,vec2}};
mod compose; use compose::BoxFn;
mod algebra; use algebra::Idx;
mod mesh; use mesh::{Mesh,Equation,Field,Operators};

struct Quantity<T=f32,const M:Mesh> { S : Field<T,M>, A : [Field<T,M>; 2] }
impl<T:Zero, const M:Mesh> Quantity<T,M> { fn new<S0:Fn(uint2)->T>(s0: S0) -> Self { Self{S:algebra::collect(|i|s0(mesh::mesh::<M>(i))), A:Zero::zero() } } }
impl<T,const M:Mesh> Operators<M> for Quantity<T,M> {}
fn mul(a:f32, b:f32) -> f32 { a*b }
impl<T:Zero+Copy+std::ops::Add<Output=T>+std::iter::Sum+'static, const M:Mesh> Quantity<T,M> where f32:std::ops::Mul<T,Output=T> {
    fn step<G:Fn(Idx)->T>(&mut self, E:&Equation<T,M>, δt:f32, g:&G) where T:std::ops::Sub<Output=T> {
        // 2-step Adams-Bashforth:  y[2] = y[1] + 3/2·δt·A(t[1],y[1]) - 1/2·δt·A(t[0],y[0])
        let b : Field<T,M> = {
            let BS = ((E.B)())(&self.S);
            algebra::collect(|i| BS(i) + (mul(3.,δt)/2.)*self.A[1][i] - (δt/2.)*self.A[0][i] + g(i))
        };
        self.S = E.A.solve(b)
    }
    fn advect(&mut self, U:&xy<Field<f32,M>>) where Self:Operators<M> {
        self.A.swap(0, 1);
        //self.A[1] = algebra::collect((&U.x*&Self::operator::<_,{Self::Dx}> + &U.y*&Self::operator::<_,{Self::Dy}>)(&self.S));
        self.A[1] = algebra::collect((&U.x*&Self::operator(Self::Dx) + &U.y*&Self::operator(Self::Dy))(&self.S));
    }
}

mod BoundaryCondition {
    use framework::core::{abs,mask};
    pub fn constant(m:u32,p:u32, d:i32) -> f32 { assert!(p==0||p==m-1); mask(d==0, 1.) }
    fn kernel<const C:[f32;3], const SYM:f32>(m:u32,p:u32,d:i32) -> f32 {
        if (-(C.len() as i32)+1..C.len() as i32).contains(&d) {
            if p==0 { C[d as usize] }
            else if p==m-1 { SYM*C[(m-1-abs(d) as u32) as usize] }
            else { panic!(); }
        } else { 0. }
    }
    pub fn derivative(m:u32,p:u32, d:i32) -> f32 { kernel::<{[-3.,4.,-1f32]},-1f32>(m,p,d)/2. } // Boundary condition on derivative
    pub fn Thom(m:u32,p:u32, d:i32) -> f32 { kernel::<{[0.,-8.,1f32]},1f32>(m,p,d) } // Thom boundary condition
}
use BoundaryCondition::*;

struct System<const M:Mesh> {
    T : Equation<f32,M>, // Temperature
    ω : Equation<f32,M>, // Vorticity
        ωT : f32, // Boussinesq approximation in buoyancy-driven flows
    φ : Equation<f32,M>, // Stream function (u=∇×φ)
    C : Equation<vec2,M>, // Color (visualization)
}

impl<const M:Mesh> Operators<M> for System<M> {}
impl<const M:Mesh> System<M> where Self:Operators<M> {
    fn new(δt : f32, Pr : f32, Ra : f32) -> Self {
        //let I=||BoxFn::new(Self::I); let P=||BoxFn::new(Self::P); let Δ=||BoxFn::new(Self::Δ);
        let I=Self::I; let P=Self::P; let Δ=Self::Δ;
        Self{
            T: Self::eq(P()      - (δt/2.)*Δ() + Self::BC::<constant,derivative>,move||box(      Self::P() + (δt/2.)*Self::Δ()     )),
            ω: Self::eq(I()  - (Pr*δt/2.)*Δ()                                                     ,move||box(      Self::P() + (Pr*δt/2.)*Self::Δ())), ωT: Ra*Pr/2.,
            φ:  Self::eq(I()-P()             + Δ()                                                  ,move||box(-1.*Self::P()                         )),
            C: Self::eq(P() - (Pr*δt/2.)*Δ() + Self::BC::<constant,constant>,move||box(      Self::P() + (Pr*δt/2.)*Self::Δ())),
        }
    }
}

struct State<const M:Mesh> {
    φ: [Field<f32,M>; 2],
    T : Quantity<f32,M>,
    ω : Quantity<f32,M>,
    C : Quantity<vec2,M>
}
impl<const M:Mesh> State<M> where System<M>:Operators<M> {
    fn C0(p:uint2) -> vec2 { p.as_f32() / (M - 1.into()).as_f32() }
    fn new() -> Self { Self{φ:Zero::zero(), T:Quantity::new(|_|0.), ω:Quantity::new(|_|0.), C:Quantity::new(Self::C0) } }
    fn C_BC(i:Idx) -> vec2 { let p = mesh::mesh::<M>(i); mask(System::border(p), Self::C0(p)) }
    fn update(&mut self, system : &System<M>, δt : f32) where Quantity<f32,M>:Operators<M>, Quantity<vec2,M>:Operators<M> {
        // Solves implicit evolution
        let Self{φ, T, ω, C} = self;
        let T0 : Field<f32,M> = T.S.clone(); // 2T½ = T0 + T1
        T.step(&system.T, δt, &|_|0.);
        // Ra·Pr·∂x/2·2T
        ω.step(&system.ω, δt, &(system.ωT*(System::operator(System::Dx))(&|i|T0[i]+T.S[i]) + System::operator(System::BC::<Thom,Thom>)(&|i|2.*φ[1][i]-φ[0][i])));
        C.step(&system.C, δt, &Self::C_BC);
        // Evaluates explicit advection
        φ.swap(0, 1); φ[1]=system.φ.A.solve(((system.φ.B)())(&ω.S));
        let U = &xy{x:algebra::collect(System::operator(System::Dy)(&φ[1])), y: algebra::collect(-System::operator(System::Dx)(&φ[1]))}; // cross(D)*(φ,φ)
        T.advect(U); ω.advect(U); C.advect(U);
   }
}

struct Parameters {Pr : f32, Ra : f32}
mod SI;
fn parameters() -> Parameters {
    use crate::SI::*;
    let η : DynamicViscosity = 2e-5 |Pa·s; //kg/m/s
    let ρ : MassDensity = 1.2 |kg_m3; //kg/m³
    let α : ThermalDiffusivity = 2e-5 |m2_s; //m²/s
    let ΔT : Temperature = 1.0 |K; // Difference between walls
    let β : ThermalExpansion = 1./(300.|K); //K¯¹: Ideal gas
    let g : Acceleration = 9.8 |m_s2; //m/s²: Gravity
    let L : Length = 0.1 |m; // Box length
    //let t : Time = sq(L)/α; //L²/α s: Box thermal diffusion time
    let ν : KinematicViscosity = η/ρ; //m²/s
    Parameters{Pr: (ν/α).into(), //Prandtl: kinematic/thermal diffusivity
                        Ra: ((ΔT*β*g*cb(L))/(ν*α)).into() //ΔT·β·g·L³/(ν·α) Rayleigh: Laminar/turbulent flow
    }
}

fn main() { //where System<M>:Operators<M> {
    const M : Mesh = Mesh{x:128, y:128};
    let Parameters{Pr:_,Ra} = parameters();
    let _δt : f32 = 1./(max(M.x,M.y) as f32*sqrt(Ra)); //R·√Ra 1: Temporal resolution
    //let system = System::new(δt, Pr, Ra);
    //let mut state = State::new();
    //state.update(&system, δt);
    /* subplot(position, size, 4, 0, Cw, Mx, My, "Vorticity ω"_);
        subplot(position, size, 4, 1, Ux, Uy, Mx, My, "Velocity u"_);
        subplot(position, size, 4, 2, positiveToImage(Ct, Mx, My), "Temperature T"_);
        subplot(position, size, 4, 3, toImage(CTx, CTy, Mx, My), "Advection"_);*/
}
