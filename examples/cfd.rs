#![allow(incomplete_features,non_snake_case)]#![feature(const_generics,non_ascii_idents,box_syntax)]use simulation::*;

struct Quantity<T=f32,const M:Mesh> { u : Field<T,M>, f : [Field<T,M>; 2] }
impl<T:Zero, const M:Mesh> Quantity<T,M> { fn new<U0:Fn(uint2)->T>(u0: U0) -> Self { Self{u:mesh::collect(u0), f:Zero::zero() } } }
impl<T:Copy+Add+Sum, const M:Mesh> Quantity<T,M> where f32:Mul<T> {
    fn step<G:Fn(Idx)->T>(&mut self, equation:&Equation<T,M>, δt:f32, g:&G) where T:Sub {
        // 2-step Adams-Bashforth:  y[2] = y[1] + 3/2·δt·f(t[1],y[1]) - 1/2·δt·f(t[0],y[0])
        let b : Field<T,M> = {
            fn mul(a:f32, b:f32) -> f32 { a*b }
            collect(|i| ((equation.B)())(&self.u)(i) + (mul(3.,δt)/2.)*self.f[1][i] - (δt/2.)*self.f[0][i] + g(i))
        };
        self.u = equation.A.solve(b)
    }
    fn advect(&mut self, v : &xy<Field<f32,M>>) {
        self.f.swap(0, 1);
        self.f[1] = collect((&v.x*&Dx!() + &v.y*&Dy!())(&self.u));
    }
}

use boundary_condition::*;
mod boundary_condition {
    use super::*;
    pub fn constant<const M:u32>(p:u32, d:i32) -> f32 { assert!(p==0||p==M-1); mask(d==0, 1.) }
    fn kernel<const C:[f32;3], const SYM:f32>(m:u32,p:u32,d:i32) -> f32 {
        if (-(C.len() as i32)+1..C.len() as i32).contains(&d) {
            if p==0 { C[d as usize] }
            else if p==m-1 { SYM*C[(m-1-abs(d) as u32) as usize] }
            else { panic!(); }
        } else { 0. }
    }
    pub fn derivative<const M:u32>(p:u32, d:i32) -> f32 { kernel::<{[-3.,4.,-1f32]},-1f32>(M,p,d)/2. } // Boundary condition on derivative
    pub fn Thom<const M:u32>(p:u32, d:i32) -> f32 { kernel::<{[0.,-8.,1f32]},1f32>(M,p,d) } // Thom boundary condition
}

struct System<const M:Mesh> {
    δt : f32,
    T : Equation<f32,M>, // Temperature
    ω : Equation<f32,M>, // Vorticity
        ωT : f32, // Boussinesq approximation in buoyancy-driven flows
    φ : Equation<f32,M>, // Stream function (u=∇×φ)
    C : Equation<vec2,M>, // Color (visualization)
}

impl<const M:Mesh> System<M> {
    fn new(δt : f32, Pr : f32, Ra : f32) -> Self {
        Self{δt,
            T: Equation::new(P!() - (    δt/2.)*Δ!() + Matrix::new(BC!(constant,derivative)),move||box(      P!() + (     δt/2.)*Δ!()     )),
            ω: Equation::new(I!() - (Pr*δt/2.)*Δ!()                                                               ,move||box(      P!() + (Pr*δt/2.)*Δ!())), ωT: Ra*Pr/2.,
            φ: Equation::new(I!()-P!()          + Δ!()                                                              ,move||box(-1.*P!()                         )),
            C: Equation::new(P!() - (Pr*δt/2.)*Δ!() + Matrix::new(BC!(constant,constant)),move||box(      P!() + (Pr*δt/2.)*Δ!())),
        }
    }
}

struct State<const M:Mesh> {
    φ: [Field<f32,M>; 2],
    T : Quantity<f32,M>,
    ω : Quantity<f32,M>,
    C : Quantity<vec2,M>
}

impl<const M:Mesh> State<M> {
    fn C0(p:uint2) -> vec2 { p.as_f32() / (M - 1.into()).as_f32() }
    fn new() -> Self { Self{φ:Zero::zero(), T:Quantity::new(|_|0.), ω:Quantity::new(|_|0.), C:Quantity::new(Self::C0) } }
    fn C_BC(i:Idx) -> vec2 { let p = mesh::<M>(i); mask(border::<M>(p), Self::C0(p)) }
    fn step(&mut self, system : &System<M>) {
        // Solves implicit evolution
        let Self{φ, T, ω, C} = self;
        let T0 : Field<f32,M> = T.u.clone(); // 2T½ = T0 + T1
        T.step(&system.T, system.δt, &|_|0.);
        ω.step(&system.ω, system.δt, &(system.ωT*Dx!()(&|i|T0[i]+T.u[i]) + operator_BC!(Thom,Thom)(&|i|2.*φ[1][i]-φ[0][i])));
        C.step(&system.C, system.δt, &Self::C_BC);
        // Evaluates explicit advection
        φ.swap(0, 1); φ[1]=system.φ.A.solve(((system.φ.B)())(&ω.u));
        let v = &xy{x:collect(Dy!()(&φ[1])), y: collect(-Dx!()(&φ[1]))};
        T.advect(v); ω.advect(v); C.advect(v);
   }
}

struct Simulation<const M:Mesh> {
    system : System<M>,
    state : State<M>
}
impl<const M:Mesh> Simulation<M> { fn new(δt : f32, Pr : f32, Ra : f32) -> Self { Self{system: System::new(δt, Pr, Ra), state: State::new()} } }

impl<const M:Mesh> Solution for Simulation<M> {
    fn current(&self) -> view::Field { view::Field{size: M, data: &self.state.T.u.0} }
    fn step(&mut self) { self.state.step(&self.system); }
}

fn main() -> Result {
    const M:Mesh = xy{x:2,y:2};

    struct Parameters {Pr : f32, Ra : f32}
    fn parameters() -> Parameters {
        use SI::*;
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
    let Parameters{Pr,Ra} = parameters();
    let δt : f32 = 1./(max(M.x,M.y) as f32*sqrt(Ra)); //R·√Ra 1: Temporal resolution
    run(View!(Simulation::<M>::new(δt, Pr, Ra)))
}
