use framework::{core::{Result,abs}, vector::{xy,uint2,size2}, image::{IntoPixelIterator, bgra8}, window::{Widget,Target}};

pub struct Field<'a> {pub size:size2, pub data:&'a [f32] }
impl std::ops::Index<uint2> for Field<'_> { type Output=f32; fn index(&self, xy{x,y}:uint2) -> &Self::Output { &self.data[(y*self.size.x+x) as usize] } }

pub trait Solution {
    fn current(&self) -> Field<'_>;
    fn step(&mut self);
}

pub struct View<T>(pub T);
impl<T:Solution> Widget for View<T> {
    fn size(&mut self, size : size2) -> size2 { size }
    fn render(&mut self, target : &mut Target) -> Result {
        let field = self.0.current();
        for (p, pixel) in target.pixels() {
            let c = field[p*(field.size-1.into())/(target.size-1.into())];
             //assert!(c >= 0. && c<= 1., c);
            let a = framework::image::sRGB::sRGB(f32::min(abs(c),1.));
            *pixel = if c>0. { bgra8{b : 0, g : a, r : a, a : 0xFF} } else { bgra8{b : a, g : a, r : 0, a : 0xFF} };
        }
        self.0.step();
        Ok(())
    }
}

#[macro_export] macro_rules! View { (Simulation::<M>::new($($args:expr),+)) => ({
        struct Instance<const MX:u32, const MY:u32>(Simulation::<{xy{x:MX,y:MY}}>);
        impl<const MX:u32, const MY:u32> Solution for Instance<MX,MY> {
            fn current(&self) -> $crate::view::Field<'_> { self.0.current() }
            fn step(&mut self) { self.0.step() }
        }
        View(Instance::<{M.x},{M.y}>(Simulation::new($($args),+)))
})}
