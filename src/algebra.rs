use framework::core::Zero;

trait Dot { fn dot<A:Fn(u32)->f32, B:Fn(u32)->f32>(self, a : A, b : B) -> f32; }
impl Dot for std::ops::Range<u32> { fn dot<A:Fn(u32)->f32, B:Fn(u32)->f32>(self, a : A, b : B) -> f32 { self.map(|k| a(k) * b(k)).sum() } }
impl Dot for std::ops::RangeTo<u32> { fn dot<A:Fn(u32)->f32, B:Fn(u32)->f32>(self, a : A, b : B) -> f32 { (0..self.end).dot(a,b) } }

pub trait SizedVector<T=f32,const N:u32> = Fn(u32)->T; //Fn<(u32,),Output=T>;
pub trait Vector<T=f32,const N:u32> = SizedVector<T,N>;

trait DotN { fn dot<A:Fn(u32)->f32, B:Fn(u32)->f32,const N:u32>(self, a : A, b : B) -> f32; }
impl DotN for std::ops::RangeFrom<u32> { fn dot<A:SizedVector<f32,N>, B:SizedVector<f32,N>, const N:u32>(self, a : A, b : B) -> f32 { (self.start..N).dot(a,b) } }
pub fn dot<A:SizedVector<N>, B:SizedVector<N>, const N:u32>(a : A, b : B) -> f32 { (..).dot(a,b) }
//pub fn dot<A:SizedVector<N>, B:SizedVector<N>, N>(a : A, b : B) -> f32 { (0..N).dot(a,b) }

pub struct Array<T=f32, const N:u32>([T; N as usize]);
pub type DenseVector<T=f32, const N:u32> = Array<T,N>;

impl<T, const N:u32>  std::ops::Index<u32> for Array<T,N> { type Output = T; fn index(&self, i : u32) -> &Self::Output { &self.0[i as usize] } }
impl<T, const N:u32>  std::ops::IndexMut<u32> for Array<T,N> { fn index_mut(&mut self, i : u32) -> &mut Self::Output { &mut self.0[i as usize] } }

impl<T:Copy,const N:u32> FnOnce<(u32,)> for Array<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (u32,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:u32> FnMut<(u32,)> for Array<T,N> { extern "rust-call" fn call_mut(&mut self, args: (u32,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:u32> Fn<(u32,)> for Array<T,N> { extern "rust-call" fn call(&self, args: (u32,)) -> Self::Output { self[args.0] } }

pub fn collect<T, F:Vector<T,N>, const N:u32>(f : F) -> Array<T, N> { Array(framework::core::array::collect(|i|f(i as u32))) }
impl<T:Zero, const N : u32> Zero for Array<T,N> { fn zero() -> Self { collect(|_|Zero::zero()) } }

pub trait Matrix<const N:u32> = Fn(u32,u32)->f32;

impl<T,const N:u32> std::ops::Mul<DenseVector<T,N>> for dyn Matrix<N> {  
    type Output = DenseVector<T,N>;
    fn mul(self, b:DenseVector<T,N>) -> Self::Output { collect(|i| dot(self.line(i), b)) }
}
use std::ops::Mul;
//impl<T,const N:u32> Mul<DenseVector<T,N>> for Box<dyn Mul<DenseVector<T,N>>> {  
//impl<T,const N:u32> Mul<DenseVector<T,N>> for Box<dyn Mul<DenseVector<T,N>,Output=DenseVector<T,N>>> {  
impl<T,const N:u32> Mul<DenseVector<T,N>> for Box<dyn Matrix<N>> {  
    type Output = <dyn Matrix<N> as std::ops::Mul<DenseVector<T,N>>>::Output;
    fn mul(self, b:DenseVector<T,N>) -> Self::Output { self.mul(b) }
}

pub struct SparseVector<T=f32,const N : u32> {
    indices : Vec<u32>,
    values : Vec<T>
}
impl<T,const N : u32> Default for SparseVector<T,N> { fn default() -> Self { Self{indices:Default::default(),values:Default::default()} } } // derive[Default] req T:Default for no reason

impl<T:Copy,const N : u32> SparseVector<T,N> {
    fn get(&self, i : u32) -> Option<T> { Some( self.values[ self.indices.binary_search(&i).ok()? ] ) }
    fn set(&mut self, i : u32, v : T) {
        framework::assert!(!self.indices.contains(&i), i, &self.indices);
        self.indices.push(i); self.values.push(v);
    }
}

impl<T:Copy+Zero,const N:u32> FnOnce<(u32,)> for SparseVector<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (u32,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:u32> FnMut<(u32,)> for SparseVector<T,N> { extern "rust-call" fn call_mut(&mut self, args: (u32,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:u32> Fn<(u32,)> for SparseVector<T,N> { extern "rust-call" fn call(&self, args: (u32,)) -> Self::Output { self.get(args.0).unwrap_or(T::zero()) } }

pub struct CSC<T=f32,const N : u32> { columns : [SparseVector<T,N>; N as usize] }
impl<T,const N : u32> Zero for CSC<T,N> { fn zero() -> Self { Self{columns: framework::core::array::collect(|_|Default::default())} } }

impl<T:Copy+Zero,const N : u32> CSC<T,N> {
    fn column(&self, j : u32) -> &SparseVector<T,N> { &self.columns[j as usize] }
    fn line(&self, i : u32) -> impl Vector<T,N> + '_ { move |j| self(i, j) }
    fn get(&self, i : u32, j : u32) -> Option<T> { self.column(j).get(i) }
    fn set(&mut self, i : u32, j : u32, v : T) { framework::assert!(i < N && j < N, i, j); self.columns[j as usize].set(i, v); }
}

impl<T:Copy+Zero,const N:u32> FnOnce<(u32,u32)> for CSC<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (u32,u32)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:u32> FnMut<(u32,u32)> for CSC<T,N> { extern "rust-call" fn call_mut(&mut self, args: (u32,u32)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:u32> Fn<(u32,u32)> for CSC<T,N> { extern "rust-call" fn call(&self, args: (u32,u32)) -> Self::Output { self.get(args.0, args.1).unwrap_or(T::zero()) } }

pub type LU<const N:u32> = CSC<f32,N>;
impl<const N:u32> LU<N> {
    pub fn new<M:Matrix<N>>(A : M) -> Self {
        let mut LU = CSC::zero();
        for i in 0..N {
            for j in i..N { LU.set(i, j, A(i, j) - (..i).dot(LU.line(i), LU.column(j))); }
            for j in i+1..N { LU.set(j, i,  (1. / LU(i, i)) * (A(j, i) - (..i).dot(LU.line(i), LU.column(j)))); }
        }
        LU // L+U-I
    }
    pub fn solve<B:SizedVector<f32,N>>(&self, b:B) -> Array<f32,N> {
        let mut y : Array<f32,N> = Zero::zero(); // Ly = b
        for i in 0..N { y[i] = b(i) - (..i).dot(self.line(i), &y); }
        let mut x : Array<f32,N> = Zero::zero(); // Ux = y
        for i in N-1..=0 { x[i] = (1. / self(i, i)) * (y[i] - (i+1..N).dot(self.line(i), &x)); }
        x
    }
}
