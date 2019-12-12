#![allow(incomplete_features)]#![feature(const_generics,non_ascii_idents)]#![allow(non_snake_case)]

macro_rules! letref { ( $self:expr => $($field:ident)+ ) => { $( let ref $field = $self.$field; )+ } }
macro_rules! letmut { ( $self:expr => $($field:ident)+ ) => { $( let ref mut $field = $self.$field; )+ } }

#[derive(Default)]
pub struct Vector<const N : u32>();
impl<const N : u32> Vector<N> {
    //fn get(&self, i : u32) -> f32 { unimplemented!() }
    fn set(&mut self, _i : u32, _v : f32) -> f32 { unimplemented!() }
}

impl<const N : u32> std::ops::Mul<&Vector<N>> for f32 {
    type Output = Vector<N>;
    fn mul(self, _: &Vector<N>) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Mul<Vector<N>> for f32 { type Output = Vector<N>; fn mul(self, b: Vector<N>) -> Self::Output { self*&b } }

impl<const N : u32> std::ops::Add<&Vector<N>> for &Vector<N> {
    type Output = Vector<N>;
    fn add(self, _: &Vector<N>) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Add<&Vector<N>> for Vector<N> { type Output = Vector<N>; fn add(self, b: &Vector<N>) -> Self::Output { &self+b } }
//impl<const N : u32> std::ops::Add<Vector<N>> for &Vector<N> { type Output = Vector<N>; fn add(self, b: Vector<N>) -> Self::Output { &b+self } }
impl<const N : u32> std::ops::Add<Vector<N>> for Vector<N> { type Output = Vector<N>; fn add(self, b: Vector<N>) -> Self::Output { self+&b } }

impl<const N : u32> std::ops::Sub<&Vector<N>> for Vector<N> {
    type Output = Vector<N>;
    fn sub(self, _: &Vector<N>) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Sub<Vector<N>> for Vector<N> { type Output = Vector<N>; fn sub(self, b: Vector<N>) -> Self::Output { self-&b } }

impl<const N : u32> std::ops::Mul<Vector<N>> for &Vector<N> {
    type Output = Vector<N>;
    fn mul(self, _: Vector<N>) -> Self::Output { unimplemented!() }
}

#[derive(Default)]
pub struct Matrix<const N : u32>();
impl<const N : u32> Matrix<N> {
    //fn get(&self, i : u32, j : u32) -> f32 { unimplemented!() }
    fn set(&mut self, _i : u32, _j : u32, _v : f32) -> f32 { unimplemented!() }

    fn identity() -> Self { unimplemented!() }
}
macro_rules! set { ($self:ident $base:expr => $( $M:ident [ $($i:expr),+ ] = $v:expr ; )* ) => { let base=$base as i32; $( $self.$M.set($((base+$i as i32)as u32),+,$v); )* } }

// fixme: expression template metaprogramming would be more runtime efficient
impl<const N : u32> std::ops::Add<&Matrix<N>> for &Matrix<N> {
    type Output = Matrix<N>;
    fn add(self, _: &Matrix<N>) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Add<&Matrix<N>> for Matrix<N> { type Output = Matrix<N>; fn add(self, B: &Self) -> Self::Output { &self+B } } // fixme: inplace
impl<const N : u32> std::ops::Add<Matrix<N>> for &Matrix<N> { type Output = Matrix<N>; fn add(self, B: Matrix<N>) -> Self::Output { B+self } }

impl<const N : u32> std::ops::Sub for &Matrix<N> {
    type Output = Matrix<N>;
    fn sub(self, _: Self) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Sub<&Matrix<N>> for Matrix<N> { type Output = Matrix<N>; fn sub(self, B: &Matrix<N>) -> Self::Output { &self-B } }
impl<const N : u32> std::ops::Sub<Matrix<N>> for &Matrix<N> { type Output = Matrix<N>; fn sub(self, B: Matrix<N>) -> Self::Output { self-&B } }
impl<const N : u32> std::ops::Sub for Matrix<N> { type Output = Matrix<N>; fn sub(self, B: Matrix<N>) -> Self::Output { self-&B } }

impl<const N : u32> std::ops::Mul<&Matrix<N>> for f32 {
    type Output = Matrix<N>;
    fn mul(self, _: &Matrix<N>) -> Self::Output { unimplemented!() }
}

impl<const N : u32> std::ops::Mul<&Vector<N>> for &Matrix<N> {
    type Output = Vector<N>;
    fn mul(self, _: &Vector<N>) -> Self::Output { unimplemented!() }
}
impl<const N : u32> std::ops::Mul<Vector<N>> for &Matrix<N> { type Output = Vector<N>; fn mul(self, b: Vector<N>) -> Self::Output { self*&b } }

mod umfpack {
    #![allow(dead_code,non_camel_case_types,non_upper_case_globals,improper_ctypes)]
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
    #[derive(Default)] pub struct UMFPACK<const N : u32>();
    pub fn factorize<const N : u32>(_A : crate::Matrix<N>) -> UMFPACK<N> { unimplemented!(); }
    impl<const N : u32> UMFPACK<N> {
        pub fn solve(&self, _b : crate::Vector<N>) -> crate::Vector<N> { unimplemented!(); }
    }
}

#[allow(non_snake_case)]
mod SI {
    #![allow(non_camel_case_types)]
    /// Quantities and units operators

    pub trait F32 {
        fn unwrap(self) -> f32;
        fn wrap(value : f32) -> Self;
    }

    #[derive(Clone,Copy)] pub struct Quantity<const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32>(f32);
    impl<const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> F32 for Quantity<A0,A1,A2,A3>{
        fn unwrap(self) -> f32 { self.0 }
        fn wrap(value : f32) -> Self { Self(value) }
    }

    pub struct Unit<Q>(std::marker::PhantomData<Q>);
    pub const fn unit<Q>() -> Unit<Q> { Unit(std::marker::PhantomData) }
    impl<Q:F32> std::ops::BitOr<Unit<Q>> for f32 { type Output = Q; fn bitor(self, _: Unit<Q>) -> Self::Output { Q::wrap(self) } }

    // quantity · quantity
    pub trait Mul<Q> { type Output : F32; }
    impl<Q:F32> Mul<Quantity<0,0,0,0>> for Q { type Output = Q; }
    impl<Q:F32+NotUnitless> Mul<Q> for Quantity<0,0,0,0> { type Output = Q; }
    impl Mul<Quantity<0,1,0,0>> for Quantity<0,1,0,0> { type Output = Quantity<0,2,0,0>; }
    impl Mul<Quantity<0,1,0,0>> for Quantity<0,2,0,0> { type Output = Quantity<0,3,0,0>; }
    impl Mul<Quantity<0,0,0,-1>> for Quantity<0,0,0,1> { type Output = Quantity<0,0,0,0>; }
    impl Mul<Quantity<0,3,0,0>> for Quantity<-2,1,0,0> { type Output = Quantity<-2,4,0,0>; }
    impl Mul<Quantity<-1,2,0,0>> for Quantity<-1,2,0,0> { type Output = Quantity<-2,4,0,0>; }

    impl<B : F32, const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> std::ops::Mul<B> for Quantity<A0,A1,A2,A3> where Self:Mul<B> {
        type Output = <Self as Mul<B>>::Output;
        fn mul(self, b: B) -> Self::Output { Self::Output::wrap(self.unwrap()*b.unwrap()) }
    }

    // quantity / quantity
    pub trait Div<Q> { type Output : F32; }
    impl<Q> Div<Q> for Q { type Output = Quantity<0,0,0,0>; }
    impl<Q:F32+NotUnitless> Div<Quantity<0,0,0,0>> for Q { type Output = Q; }
    impl Div<Quantity<0,0,0,1>> for Quantity<0,0,0,0> { type Output = Quantity<0,0,0,-1>; }
    impl Div<Quantity<-1,2,0,0>> for Quantity<0,2,0,0> { type Output = Quantity<1,0,0,0>; }
    impl Div<Quantity<0,-3,1,0>> for Quantity<-1,-1,1,0> { type Output = Quantity<-1,2,0,0>; }

    impl<B:F32, const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> std::ops::Div<B> for Quantity<A0,A1,A2,A3> where Self:Div<B> {
        type Output = <Self as Div<B>>::Output;
        fn div(self, b: B) -> Self::Output { Self::Output::wrap(self.unwrap()/b.unwrap()) }
    }

    pub type Unitless = Quantity<0,0,0,0>;

    // unitless · quantity
    impl<const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> std::ops::Mul<Quantity<A0,A1,A2,A3>> for f32 where Quantity<A0,A1,A2,A3>:NotUnitless {
        type Output = Quantity<A0,A1,A2,A3>;
        fn mul(self, b: Quantity<A0,A1,A2,A3>) -> Self::Output { Unitless::wrap(self)*b }
    }

    // quantity / unitless
    impl<const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> std::ops::Div<f32> for Quantity<A0,A1,A2,A3> where Self:NotUnitless {
        type Output = Self;
        fn div(self, b: f32) -> Self { self/Unitless::wrap(b) }
    }

    // unitless / quantity
    impl<const A0 : i32, const A1 : i32, const A2 : i32, const A3 : i32> std::ops::Div<Quantity<A0,A1,A2,A3>> for f32 where Unitless:Div<Quantity<A0,A1,A2,A3>> {
        type Output = <Unitless as Div<Quantity<A0,A1,A2,A3>>>::Output;
        fn div(self, b: Quantity<A0,A1,A2,A3>) -> Self::Output { Unitless::wrap(self)/b }
    }

     // f32 · unitless
    impl std::ops::Mul<Unitless> for f32 { type Output = f32; fn mul(self, b: Unitless) -> Self::Output { self*b.unwrap() } }
    //  unitless · f32
    impl std::ops::Mul<f32> for Unitless { type Output = f32; fn mul(self, b: f32) -> Self::Output { self.unwrap()*b } }
    // unitless / f32
    impl std::ops::Div<f32> for Unitless { type Output = f32; fn div(self, b: f32) -> Self::Output { self.unwrap()/b } }

    pub trait NotUnitless {}
    macro_rules! quantity_unit { ( [ $($dimensions:expr),+ ] $unit:ident $quantity:ident  ) => {
            #[allow(non_camel_case_types)] pub type $quantity = Quantity<$($dimensions),+>;
            impl NotUnitless for $quantity {}
            #[allow(dead_code,non_upper_case_globals)] pub const $unit : Unit<$quantity> = unit();
    } }

    // s m kg K
    quantity_unit!([1,0,0,0] s Time);
    quantity_unit!([0,1,0,0] m Length );
    quantity_unit!([0,0,1,0] kg Mass);
    quantity_unit!([0,0,0,1] K Temperature);
    quantity_unit!([0, 2,0,0] m2 Area);
    quantity_unit!([0, 3,0,0] m3 Volume);
    quantity_unit!([-2,1,0,0] m_s2 Acceleration);
    quantity_unit!([-1,2,0,0] m2_s Diffusivity);
    quantity_unit!([0,-3,1,0] kg_m3 Mass_density);
    quantity_unit!([-1,-1,1,0] Pa·s Dynamic_viscosity); //kg/m/s
    quantity_unit!([0,0,0,-1] _K Thermal_expansion);

    pub type Prandtl = Unitless;
    pub type Rayleigh = Unitless;
    //pub type Timestep = Unitless;
    pub type Thermal_diffusivity = Diffusivity; // m²/s
    pub type Kinematic_viscosity = Diffusivity; // m²/s

    pub use framework::{sq, cb};
    pub fn sqrt(x: Unitless) -> Unitless { Unitless::wrap(x.0.sqrt()) }
}

#[allow(non_snake_case, non_upper_case_globals)]
mod Box {
    pub const R : u32 = 128; // Resolution
    pub const Mx : u32 = R+1;
    pub const My : u32 = R+1;
    pub const N : u32 = Mx*My; // Mesh vertex count
    pub const dx : f32 = 1./Mx as f32; // Spatial resolution
    pub const dy : f32 = 1./My as f32; // Spatial resolution
}

use crate::Box::*;
use umfpack::*;
#[derive(Default)] struct Implicit<const N : u32> { L_T : UMFPACK<N>, L_ω : UMFPACK<N>, L_φ : UMFPACK<N>, L_𝜓 : UMFPACK<N> } // Factorized left-hand operator (implicit)
#[derive(Default)] struct Explicit<const N : u32> { R_T : Matrix<N>, R_ω : Matrix<N>, R_φ : Matrix<N>, R_𝜓 : Matrix<N> } // Right-hand operators (explicit)
#[allow(non_snake_case)]#[derive(Default)] struct System<const N : u32> {
    P : Matrix<N>, Δ : Matrix<N>, dx : Matrix<N>, dy : Matrix<N>, // Elementary operators: interior points projector, Laplacian Δ and partial derivatives ∂x, ∂y
    BC_T : Matrix<N>, BC_ω : Matrix<N>, BC_𝜓 : Matrix<N>, // Boundary conditions
    G_T : Vector<N>, G_𝜓x : Vector<N>, G_𝜓y : Vector<N>, // Right-hand vectors (BC and source field)
    implicit : Implicit::<N>,
    explicit : Explicit::<N>,
}
impl<const N : u32> System<N>{
    fn new(δt : f32, Pr : f32) -> Self { let mut system = System::default(); system.initialize(δt, Pr); system }
    fn initialize(&mut self, δt : f32, Pr : f32) {
        // Elementary matrices (interior points projector, Laplacian Δ and partial derivatives ∂x, ∂y)
        for x in 1..Mx-1 { for y in 1..My-1 { set!{self y*Mx+x =>
                P[0, 0] = 1.;
                                                Δ[0, -(Mx as i32)] = 1./(dy*dy);
                Δ[0, -1] = 1./(dx*dx); Δ[0, 0] = -2.*(1./(dx*dx)+1./(dy*dy)); Δ[0, 1] = 1./(dx*dx);
                                                Δ[0, Mx] = 1./(dy*dy);
                dx[0, -1] = -1./(2.*dx);
                dx[0,  1] =   1./(2.*dx);
                dy[0, -(Mx as i32)] = -1./(2.*dy);
                dy[0, 0+Mx] = 1./(2.*dy);
        }}}

        // Boundary condition for temperature T
        for x in 1..Mx-1 { set!{self x => // Constant derivative (Neumann) on horizontal boundaries
            // Top
            BC_T[0, 0*Mx] = -3./(2.*dy);
            BC_T[0, 1*Mx] = 4./(2.*dy);
            BC_T[0, 2*Mx] = -1./(2.*dy);
            G_T[0] = 0.; // No flux
            // Bottom
            BC_T[(My-1)*Mx, (My-3)*Mx] = 1./(2.*dy);
            BC_T[(My-1)*Mx, (My-2)*Mx] = -4./(2.*dy);
            BC_T[(My-1)*Mx, (My-1)*Mx] = 3./(2.*dy);
            G_T[0] = 0.; // No flux
        }}
        for y in 0..My { set!{self y*Mx => // Constant value (Dirichlet) on vertical boundaries
            // Left
            BC_T[0, 0] = 1.;
            G_T[0] = 0.;
            // Right
            BC_T[Mx-1, Mx-1] = 1.;
            G_T[Mx-1] = 1.;
        }}

        // Thom boundary condition for vorticity ω
        for x in 1..Mx-1 { set!{self x => // Horizontal boundaries
            // Top
            BC_ω[0, 1*Mx] = -8./(2.*dx*dx);
            BC_ω[0, 2*Mx] = 1./(2.*dx*dx);
            // Bottom
            BC_ω[(My-1)*Mx, (My-3)*Mx] = 1./(2.*dx*dx);
            BC_ω[(My-1)*Mx, (My-2)*Mx] = -8./(2.*dx*dx);
        }}
        for y in 0..My { set!{self y*Mx => // Vertical boundaries
            // Left
            BC_ω[0, 1] = -8./(2.*dx*dx);
            BC_ω[0, 2] = 1./(2.*dx*dx);
            // Right
            BC_ω[Mx-1, Mx-3] = 1./(2.*dx*dx);
            BC_ω[Mx-1, Mx-2] = -8./(2.*dx*dx);
        }}

        // Constant boundary condition for advection
        for x in 1..Mx-1 { set!{self x => // Horizontal boundaries
            // Top
            BC_𝜓[0*Mx, 0*Mx] = 1.;
            G_𝜓x[0*Mx] = x as f32/Mx as f32; G_𝜓y[0*Mx] = 0./My as f32;
            // Bottom
            BC_𝜓[(My-1)*Mx, (My-1)*Mx] = 1.;
            G_𝜓x[(My-1)*Mx] = x as f32/Mx as f32; G_𝜓y[(My-1)*Mx] = (My-1) as f32/My as f32;
        }}
        for y in 0..My { set!{self y*Mx => // Vertical boundaries
            // Left
            BC_𝜓[0, 0] = 1.;
            G_𝜓x[0] = 0./Mx as f32; G_𝜓y[0] = y as f32/My as f32;
            // Right
            BC_𝜓[Mx-1,-1] = 1.;
            G_𝜓x[Mx-1] = (Mx-1) as f32/Mx as f32; G_𝜓y[Mx-1] = y as f32/My as f32;
        }}

        let I = &Matrix::<N>::identity();
        letref!{self => P Δ BC_T BC_𝜓}
        self.implicit = Implicit::<N>{ // Left-hand side (implicit) of temperature T, vorticity ω, stream function φ and advection 𝜓 evolution equations
            L_T: factorize( BC_T + P - (δt/2.)*Δ ),
            L_ω: factorize( I - (Pr*δt/2.)*Δ ),
            L_φ: factorize( (I-P) + Δ ),
            L_𝜓: factorize( BC_𝜓 + P - (1./R as f32*δt/2.)*Δ )
        };
        self.explicit = Explicit::<N>{ // Right-hand side (explicit)
            R_T: P + (δt/2.)*Δ,
            R_ω: P + (Pr*δt/2.)*Δ,
            R_φ: -1.*P,
            R_𝜓: P + (1./R as f32*δt/2.)*Δ
        };
    }
}

#[derive(Default)]
struct State<const N: u32> {
    C_T : Vector<N>, // Current temperature T field
    C_ω : Vector<N>, // Current vorticity ω field (ω=∇×u)
    P_φ : Vector<N>, C_φ : Vector<N>, // Previous and current stream function φ field (u=∇×φ)
    C_𝜓x : Vector<N>, C_𝜓y : Vector<N>, // Current advectionfields

    PA_T : Vector<N>, CA_T : Vector<N>, // Previous and current non-linear advection term for temperature T
    PA_ω : Vector<N>, CA_ω : Vector<N>, // Previous and current non-linear advection term for vorticity ω
    PA_𝜓x : Vector<N>, CA_𝜓x : Vector<N>, // Previous and current non-linear advection
    PA_𝜓y : Vector<N>, CA_𝜓y : Vector<N>, // Previous and current non-linear advection
}
impl<const N : u32> State<N> {
    fn new() -> Self { let mut state = State::default(); state.initialize(); state }
    fn initialize(&mut self) {
        // Sets initial positions for advection
        for x in 0..Mx { for y in 0..My { set!{self y*Mx+x =>
            C_𝜓x[0] = x as f32/Mx as f32;
            C_𝜓y[0] = y as f32/My as f32;
        } } }
    }
    fn update(&mut self, system : &System<N>, δt : f32, Ra : f32, Pr : f32) {
        // Solves evolution equations
        letref!{system.implicit => L_T L_ω  L_φ  L_𝜓}
        letref!{system.explicit => R_T R_ω R_φ R_𝜓}
        letref!{system => G_T BC_ω G_𝜓x G_𝜓y}
        letref!{self => C_T C_ω C_φ C_𝜓x C_𝜓y P_φ}
        letref!{self => CA_T CA_ω CA_𝜓x CA_𝜓y}
        letref!{self => PA_T PA_ω PA_𝜓x PA_𝜓y}
        let N_T = L_T.solve(R_T*C_T   + (3.*δt/2.)*CA_T  - (δt/2.)*PA_T  + G_T);
        let N_ω = L_ω.solve(R_ω*C_ω + (3.*δt/2.)*CA_ω - (δt/2.)*PA_ω + (Ra*Pr*δt/2.)*(dx*(C_T+&N_T)) + BC_ω*(2.*C_φ-P_φ));
        let N_φ = L_φ.solve(R_φ*&N_ω);
        let N_𝜓x = L_𝜓.solve(R_𝜓*C_𝜓x   + (3.*δt/2.)*CA_𝜓x  - (δt/2.)*PA_𝜓x  + G_𝜓x);
        let N_𝜓y = L_𝜓.solve(R_𝜓*C_𝜓y   + (3.*δt/2.)*CA_𝜓y  - (δt/2.)*PA_𝜓y  + G_𝜓y);

        // Update references
        letmut!{self => C_T C_ω C_𝜓x C_𝜓y P_φ C_φ}
        *C_T=N_T; *C_ω=N_ω; *C_𝜓x=N_𝜓x; *C_𝜓y=N_𝜓y;
        use std::mem::swap;
        swap(P_φ, C_φ);
        letmut!{self => C_φ}; *C_φ=N_φ;
        letmut!{self => CA_T CA_ω CA_𝜓x CA_𝜓y}
        letmut!{self => PA_T PA_ω PA_𝜓x PA_𝜓y}
        swap(PA_T, CA_T); swap(PA_ω, CA_ω); swap(PA_𝜓x, CA_𝜓x); swap(PA_𝜓y, CA_𝜓y);

        // Computes advection for next step
        letref!{self => C_φ C_T C_ω C_𝜓x C_𝜓y}
        let ref Ux =       dy*C_φ;
        let ref Uy = -1.*dx*C_φ;
        letmut!{self => CA_T CA_ω CA_𝜓x CA_𝜓y}
        *CA_T     = Ux*(dx*C_T   )+Uy*(dy*C_T);
        *CA_ω    = Ux*(dx*C_ω  )+Uy*(dy*C_ω);
        *CA_𝜓x  = Ux*(dx*C_𝜓x)+Uy*(dy*C_𝜓x);
        *CA_𝜓y  = Ux*(dx*C_𝜓y)+Uy*(dy*C_𝜓y);
    }
}

fn main() {
    use crate::SI::*;
    let L : Length = 0.1 |m; // Box side
    let α : Thermal_diffusivity = 2e-5 |m2_s; //m²/s
    let η : Dynamic_viscosity = 2e-5 |Pa·s; //kg/m/s
    let ρ : Mass_density = 1.2 |kg_m3; //kg/m³
    let ΔT : Temperature = 1.0 |K; // Difference between walls
    let β : Thermal_expansion = 1./(300.|K); //K¯¹: Ideal gas
    let g : Acceleration = 9.8 |m_s2; //m/s²: Gravity

    //let t : Time = sq(L)/α; //L²/α s: Box thermal diffusion time
    let ν : Kinematic_viscosity = η/ρ; //m²/s
    let Pr : Prandtl = ν/α; //1: momentum/thermal diffusivity
    let Ra : Rayleigh = (ΔT*β*g*cb(L))/(ν*α); //ΔT·β·g·L³/(ν·α) 1: Laminar vs turbulent flow
    let δt : f32 = 1./(R as f32*sqrt(Ra)); //R·√Ra 1: Temporal resolution

    let system = System::<N>::new(δt, 1.*Pr);
    let mut state = State::<N>::new();
    state.update(&system, δt, 1.*Ra, 1.*Pr);
    /* subplot(position, size, 4, 0, Cw, Mx, My, "Vorticity ω"_);
        subplot(position, size, 4, 1, Ux, Uy, Mx, My, "Velocity u"_);
        subplot(position, size, 4, 2, positiveToImage(Ct, Mx, My), "Temperature T"_);
        subplot(position, size, 4, 3, toImage(CTx, CTy, Mx, My), "Advection"_);*/
}
