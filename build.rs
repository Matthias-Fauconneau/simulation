fn main() {
    println!("cargo:rustc-link-lib=umfpack");
    println!("cargo:rerun-if-changed=wrapper.h");
    bindgen::Builder::default().header("wrapper.h").parse_callbacks(Box::new(bindgen::CargoCallbacks)).generate().unwrap()
        .write_to_file(std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("bindings.rs")).unwrap();
}
