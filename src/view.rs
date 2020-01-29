use framework::{core::{Result,abs}, vector::{xy,uint2,size2}, image::{IntoPixelIterator, bgra8}, window::{Widget,Target}};

pub struct Field<'a> { size:size2, data:&'a [f32] }
//impl std::ops::Deref for Field { type Target = &[f32]; fn deref(&self) -> &Self::Target { &self.data } }
//impl<T,const M:Mesh> std::ops::Deref for Field<T,M> { type Target = algebra::Array<T,{N(M)}>; fn deref(&self) -> &Self::Target { &self.0 } }
//impl<T,const M:Mesh> std::ops::DerefMut for Field<T,M> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 } }
//pub fn collect<T, F:Fn(Idx)->T, const M:Mesh>(f : F) ->Field<T,M> { Field(algebra::collect(f)) }
//impl<T:Zero, const M:Mesh> Zero for Field<T,M> { fn zero() -> Self { Field(Zero::zero()) } }
impl std::ops::Index<uint2> for Field<'_> { type Output=f32; fn index(&self, xy{x,y}:uint2) -> &Self::Output { &self.data[(y*self.size.x+x) as usize] } }
/*impl<'a> FnOnce<(uint2,)> for Field { type Output=&'a T;  extern "rust-call" fn call_once(self, args: (uint2,)) -> Self::Output { self.call(args) } }
impl<'a> FnMut<(uint2,,)> for Field { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<'a> Fn<(uint2,,)> for Field { extern "rust-call" fn call(&self, args: (uint2,)) -> Self::Output { self[args.0] } }*/

pub trait Solution {
    fn current(&self) -> Field<'_>;
    fn step(&mut self);
}

impl Widget for dyn Solution {
    fn size(&mut self, size : size2) -> size2 { size }
    fn render(&mut self, target : &mut Target) -> Result {
        let field = self.current();
        for (p, pixel) in target.pixels() {
            let c = field[p*(field.size-1.into())/(target.size-1.into())];
             //assert!(c >= 0. && c<= 1., c);
            let a = framework::image::sRGB::sRGB(f32::min(abs(c),1.));
            *pixel = if c>0. { bgra8{b : 0, g : a, r : a, a : 0xFF} } else { bgra8{b : a, g : a, r : 0, a : 0xFF} };
        }
        self.step();
        Ok(())
    }
}

#[macro_export] macro_rules! run { (Simulation::<M>::new($($args:expr),+)) => (
    {
        struct Run<const MX:u32, const MY:u32>(Simulation<{xy{x:MX,y:MY}}>);
        use framework::window::Widget;
        impl<const MX:u32, const MY:u32> Widget for Run<MX,MY> {
            fn size(&mut self, size:framework::vector::size2) -> framework::vector::size2 { self.0.size(size) }
            fn render(&mut self, target:&mut framework::window::Target) -> Result { self.0.render(target) }
        }
        framework::window::run(Run::<{M.x},{M.y}>(Simulation::new($($args),+)))
    }
)}
