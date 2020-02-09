use framework::{core::abs,vector::{xy,uint2,size2}, image::bgra8, window::{Widget,Target}};
fn min(x: &[f32]) -> &f32 { x.iter().min_by(|a, b| a.partial_cmp(b).expect("NaN")).unwrap() }
fn max(x: &[f32]) -> &f32 { x.iter().max_by(|a, b| a.partial_cmp(b).expect("NaN")).unwrap() }

pub struct Field<'t> {pub size:size2, pub values:&'t [f32] }
impl std::ops::Index<uint2> for Field<'_> { type Output=f32; fn index(&self, xy{x,y}:uint2) -> &Self::Output { &self.values[(y*self.size.x+x) as usize] } }

pub trait Solution {
    fn current(&self) -> Field<'_>;
    fn step(&mut self);
}

pub struct View<T>(pub T);
impl<T:Solution> Widget for View<T> {
    fn size(&mut self, size : size2) -> size2 { size/2 }
    fn render(&mut self, target : &mut Target) {
        let field = self.0.current();
        let min = min(field.values);
        let max = max(field.values);
        //log!(min, max);
        let time = std::time::Instant::now();
        let size = target.size;
        target.set(|p| {
            let c = (field[p*(field.size-1.into())/(size-1.into())]-min)/(max-min);
            //assert!(c >= 0. && c<= 1., c);
            let a = framework::image::sRGB::sRGB(f32::min(abs(c),1.));
            if c>0. { bgra8{b : 0, g : a, r : a, a : 0xFF} } else { bgra8{b : a, g : a, r : 0, a : 0xFF} }
        });
        //log!(time.elapsed().as_millis()); let time = std::time::Instant::now();
        self.0.step(); // FIXME: async
        print!("{} ",time.elapsed().as_millis()); use std::io::{Write,stdout}; stdout().flush().unwrap();
        //Ok(())
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
