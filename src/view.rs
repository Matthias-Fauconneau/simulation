use framework::{log, core::{abs,div_rem}, vector::{xy,uint2,size2,vec2,norm,atan}, image::{sRGB::sRGB, bgra8}, window::{Widget,Target}, color::LCh};
trait From<T> { fn from(t: T) -> Self; }
impl From<vec2> for bgra8 { fn from(v: vec2) -> Self { std::convert::Into::into(LCh{L:norm(v), C:norm(v), h:atan(v)+std::f32::consts::PI}) } }
trait Into<T> { fn into(self) -> T; }
impl<T,S> Into<T> for S where T:From<Self> { fn into(self) -> T { T::from(self) } }
trait TryMax<T> : Iterator<Item=T> { fn try_max(self) -> Self::Item; }
impl<I,T:PartialOrd> TryMax<T> for I where Self:Iterator<Item=T> { fn try_max(self) -> Self::Item { self.max_by(|a, b| a.partial_cmp(b).expect("NaN")).unwrap() } }

pub enum Field<'t> { Positive(&'t [f32]), Scalar(&'t [f32]), XY(&'t [vec2]) }

pub trait Solution {
    fn output(&self) -> (size2, [Field<'_>; 4]);
    fn step(&mut self);
}

pub struct View<T>(pub T);
impl<T:Solution> Widget for View<T> {
    fn size(&mut self, size : size2) -> size2 { xy{x:size.y,y:size.y} }
    fn render(&mut self, target : &mut Target) {
        self.0.step(); // FIXME: async
        let (field_size, outputs) = self.0.output();
        for i in 0..outputs.len() {
            let mut target = { let size = target.size/2; let (div,rem) = div_rem(i as u32,2); target.slice_mut(xy{x:rem, y:div}*size, size) };
            let index = {let target_size=target.size; move |p:uint2| -> usize { let xy{x,y} = p*field_size/target_size; (y*field_size.x+x) as usize }};
            match outputs[i as usize] {
                Field::Positive(field) => {
                    let max = field.iter().try_max();
                    log!(i,"+",max);
                    target.set(|p| {
                        let s = sRGB( field[index(p)] / max );
                        bgra8{b:s, g:s, r:s, a:0xFF}
                    });
                }
                Field::Scalar(field) => {
                    let max = field.iter().map(|v|abs(*v)).try_max();
                    if max == 0. {continue;}
                    log!(i,"~",max);
                    target.set(|p| {
                        let v = field[index(p)];
                        let s = sRGB( abs(v) / max );
                        if v>0. {bgra8{b:0, g:s, r:s, a:0xFF}} else {bgra8{b:s, g:s, r:0, a:0xFF}}
                    });
                }
                Field::XY(field) => {
                    let max = field.iter().map(|v|norm(*v)).try_max();
                    if max == 0. {continue;}
                    log!(i,".",max);
                    target.set(|p| { self::Into::into((1./max)*field[index(p)]) });
                }
            }
        }
    }
}

#[macro_export] macro_rules! View { (Simulation::<M>::new($($args:expr),+)) => ({
        struct Instance<const MX:u32, const MY:u32>(Simulation::<{xy{x:MX,y:MY}}>);
        impl<const MX:u32, const MY:u32> Solution for Instance<MX,MY> {
            fn output(&self) -> (size2, [$crate::view::Field<'_>; 4]) { self.0.output() }
            fn step(&mut self) { self.0.step() }
        }
        View(Instance::<{M.x},{M.y}>(Simulation::new($($args),+)))
})}
