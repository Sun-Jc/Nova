//! Lightweight phase-timing hook for profiling `ppsnark::prove`.
//!
//! Feature-gated behind `phase-profile` so it never touches production builds.
//! `prove()` calls [`mark`] at each phase boundary; a bench harness installs a
//! callback via [`set_hook`] to record wall-clock time and (its own) RSS sample
//! at each boundary. The library side stays `unsafe`-free — any RSS syscall
//! lives in the harness's callback.
use std::sync::RwLock;

/// A phase-boundary callback: receives the label of the boundary just reached.
type Hook = Box<dyn Fn(&str) + Send + Sync>;

static HOOK: RwLock<Option<Hook>> = RwLock::new(None);

/// Install the phase-boundary callback. Passing `None` clears it.
pub fn set_hook(hook: Option<Hook>) {
  *HOOK.write().unwrap() = hook;
}

/// Signal that a phase boundary named `label` has been reached.
///
/// No-op when no hook is installed. Called from `ppsnark::prove` at each phase
/// boundary; the string label matches the phase names in the bench report.
pub fn mark(label: &str) {
  if let Some(hook) = HOOK.read().unwrap().as_ref() {
    hook(label);
  }
}
