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
quantity_unit!([0,-3,1,0] kg_m3 MassDensity);
quantity_unit!([-1,-1,1,0] Pa·s DynamicViscosity); //kg/m/s
quantity_unit!([0,0,0,-1] _K ThermalExpansion);

pub type ThermalDiffusivity = Diffusivity; // m²/s
pub type KinematicViscosity = Diffusivity; // m²/s

impl From<Unitless> for f32 { fn from(v: Unitless) -> Self { v.unwrap() } }
