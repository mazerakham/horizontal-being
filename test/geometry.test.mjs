// node test/geometry.test.mjs
// Tests the geometry block that ships inside index.html (extracted between the BEGIN/END GEOMETRY markers),
// against numerical integration of the horizontal lift and brute-force marching.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const html = readFileSync(join(here, '..', 'index.html'), 'utf8');
const m = /\/\/ BEGIN GEOMETRY([\s\S]*?)\/\/ END GEOMETRY/.exec(html);
if (!m) throw new Error('geometry block not found in index.html');
const HG = new Function(m[1] + '\nreturn HG;')();

let failures = 0;
function report(name, err, tol) {
  const ok = err <= tol;
  if (!ok) failures++;
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}: max error ${err.toExponential(2)} (tolerance ${tol.toExponential(0)})`);
}
let seed = 12345;
const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
const U = (a, b) => a + (b - a) * rnd();
const dist = (a, b) => Math.max(Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]), Math.abs(a[2] - b[2]));

// RK4 for the horizontal lift: (x, y, t)' = (a, b, 2(x b - y a)) with planar velocity (a, b) = dir(s).
function liftRK4(start, dirAt, len, steps) {
  let [x, y, t] = start; const h = len / steps;
  const f = (s, X, Y) => { const [a, b] = dirAt(s); return [a, b, 2 * (X * b - Y * a)]; };
  for (let i = 0; i < steps; i++) {
    const s = i * h;
    const k1 = f(s, x, y);
    const k2 = f(s + h / 2, x + h / 2 * k1[0], y + h / 2 * k1[1]);
    const k3 = f(s + h / 2, x + h / 2 * k2[0], y + h / 2 * k2[1]);
    const k4 = f(s + h, x + h * k3[0], y + h * k3[1]);
    x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    y += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    t += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
  }
  return [x, y, t];
}

// 1. closed-form probe vs ODE, local frame
{
  let err = 0;
  const kappas = [0, 1e-9, -1e-9, 1e-6, -3e-4, 0.004, -0.02];
  for (let n = 0; n < 400; n++) {
    const th = U(-Math.PI, Math.PI), k = n < kappas.length ? kappas[n] : U(-3, 3), ell = U(0, 30);
    const ode = liftRK4([0, 0, 0], (s) => [Math.cos(th + k * s), Math.sin(th + k * s)], ell, 4000);
    err = Math.max(err, dist(HG.probe(th, k, ell), ode));
  }
  report('probe(theta, kappa, l) vs RK4 horizontal lift from identity (400 probes, l up to 30, |kappa| up to 3)', err, 1e-7);
}

// 2. world composition pose . R_h(local) vs ODE started at the pose
{
  let err = 0;
  for (let n = 0; n < 300; n++) {
    const pose = { x: U(-20, 20), y: U(-20, 20), t: U(-100, 100), h: U(-Math.PI, Math.PI) };
    const th = U(-Math.PI, Math.PI), k = U(-2, 2), ell = U(0, 25);
    const ode = liftRK4([pose.x, pose.y, pose.t], (s) => [Math.cos(pose.h + th + k * s), Math.sin(pose.h + th + k * s)], ell, 4000);
    err = Math.max(err, dist(HG.worldOf(pose, HG.probe(th, k, ell)), ode));
  }
  report('worldOf(pose, probe) vs RK4 lift started at the pose (300 cases, |x|,|y| up to 20)', err, 1e-7);
}

// 3. pose integration with exact arc steps vs ODE with the same piecewise-constant inputs
{
  let err = 0;
  seed = 999;
  for (let n = 0; n < 40; n++) {
    const init = { x: U(-5, 5), y: U(-5, 5), t: U(-20, 20), h: U(-3, 3) };
    const segs = [];
    for (let j = 0; j < 30; j++) segs.push({ vf: U(-3, 3), vs: rnd() < 0.3 ? U(-2, 2) : 0, om: U(-2, 2), dt: U(0.005, 0.3) });
    let pose = { ...init };
    for (const s of segs) pose = HG.step(pose, s.vf, s.vs, s.om, s.dt);
    let [x, y, t, h] = [init.x, init.y, init.t, init.h];
    for (const s of segs) {
      const N = 400, dt = s.dt / N;
      const f = (X, Y, T, Hh) => { const a = Math.cos(Hh) * s.vf - Math.sin(Hh) * s.vs, b = Math.sin(Hh) * s.vf + Math.cos(Hh) * s.vs; return [a, b, 2 * (X * b - Y * a), s.om]; };
      for (let i = 0; i < N; i++) {
        const k1 = f(x, y, t, h);
        const k2 = f(x + dt / 2 * k1[0], y + dt / 2 * k1[1], t, h + dt / 2 * k1[3]);
        const k3 = f(x + dt / 2 * k2[0], y + dt / 2 * k2[1], t, h + dt / 2 * k2[3]);
        const k4 = f(x + dt * k3[0], y + dt * k3[1], t, h + dt * k3[3]);
        x += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        y += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
        t += dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
        h += dt / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
      }
    }
    err = Math.max(err, dist([pose.x, pose.y, pose.t], [x, y, t]), Math.abs(pose.h - h));
  }
  report('step() exact arcs vs RK4 of (x,y,t,h) under 30 random drive/strafe/turn segments (40 runs)', err, 1e-7);
}

// 4. group structure: associativity, inverse, rotation automorphism, left invariance of the gauge distance
{
  let err = 0;
  for (let n = 0; n < 2000; n++) {
    const a = [U(-9, 9), U(-9, 9), U(-50, 50)], b = [U(-9, 9), U(-9, 9), U(-50, 50)], c = [U(-9, 9), U(-9, 9), U(-50, 50)], h = U(-4, 4);
    err = Math.max(err,
      dist(HG.mul(HG.mul(a, b), c), HG.mul(a, HG.mul(b, c))),
      dist(HG.mul(a, HG.inv(a)), [0, 0, 0]),
      dist(HG.rot(h, HG.mul(a, b)), HG.mul(HG.rot(h, a), HG.rot(h, b))),
      Math.abs(HG.gauge(HG.mul(HG.inv(HG.mul(c, a)), HG.mul(c, b))) - HG.gauge(HG.mul(HG.inv(a), b))));
  }
  report('associativity, inverse, rotation automorphism, left-invariant gauge distance (2000 triples)', err, 1e-9);
}

// 5. horizontality of the world probe curve: alpha(gamma') = t' - 2(x y' - y x') = 0
{
  let err = 0;
  for (let n = 0; n < 300; n++) {
    const pose = { x: U(-10, 10), y: U(-10, 10), t: U(-40, 40), h: U(-3, 3) }, th = U(-3, 3), k = U(-2, 2), l = U(0.1, 20), e = 1e-5;
    const p = HG.worldOf(pose, HG.probe(th, k, l)), pp = HG.worldOf(pose, HG.probe(th, k, l + e)), pm = HG.worldOf(pose, HG.probe(th, k, l - e));
    const d = [(pp[0] - pm[0]) / (2 * e), (pp[1] - pm[1]) / (2 * e), (pp[2] - pm[2]) / (2 * e)];
    err = Math.max(err, Math.abs(d[2] - 2 * (p[0] * d[1] - p[1] * d[0])), Math.abs(Math.hypot(d[0], d[1]) - 1));
  }
  report('probe curves are horizontal and unit speed, alpha(gamma\') = 0 by central differences', err, 1e-5);
}

// 6. loops: CCW unit circle raises t by 4 pi; loopArea matches the shoelace area of the trail
{
  const full = HG.probe(0, 1, 2 * Math.PI);
  let err = Math.abs(full[2] - 4 * Math.PI) + Math.hypot(full[0], full[1]);
  const cw = HG.probe(0.7, -0.5, 4 * Math.PI);
  err = Math.max(err, Math.abs(cw[2] + 4 * Math.PI * 4));
  for (let n = 0; n < 50; n++) {
    let pose = { x: U(-5, 5), y: U(-5, 5), t: U(-5, 5), h: U(-3, 3) };
    const mark = { ...pose }, pts = [[pose.x, pose.y]];
    for (let j = 0; j < 600; j++) { pose = HG.step(pose, U(0, 3), 0, U(-3, 3), 0.01); pts.push([pose.x, pose.y]); }
    let sh = 0;
    for (let i = 0; i < pts.length; i++) { const a = pts[i], b = pts[(i + 1) % pts.length]; sh += a[0] * b[1] - a[1] * b[0]; }
    // trail sampled finely; chord closes the polygon, so compare to the exact loopArea with a sampling tolerance
    err = Math.max(err, Math.abs(sh / 2 - HG.loopArea(mark, pose)) > 5e-3 ? 1 : 0);
  }
  report('CCW circle lifts 4 pi r^2, CW lowers it; loopArea agrees with shoelace (50 random trails)', err, 1e-9);
}

// 7. eye hit test (the algorithm mirrored in the shader) vs brute-force dense march
{
  const L = 30;
  let probes = 0, hitsBrute = 0, idMismatch = 0, grazeMismatch = 0, maxEll = 0;
  const bad = [];
  for (let scene = 0; scene < 60; scene++) {
    const objs = [];
    const nObj = 6;
    for (let i = 0; i < nObj; i++) {
      const kind = rnd();
      let x, y, t;
      if (kind < 0.35) { x = U(-12, 12); y = U(-12, 12); t = U(-2, 2); }             // near the plane
      else if (kind < 0.6) { x = U(-2, 2); y = U(-2, 2); t = U(-60, 60); }           // above or below
      else { x = U(-15, 15); y = U(-15, 15); t = U(-120, 120); }
      objs.push({ x, y, t, r: U(0.5, 2), active: true });
    }
    for (let n = 0; n < 150; n++) {
      probes++;
      const th = U(-Math.PI, Math.PI), v = U(-1, 1), k = 2 * Math.sign(v) * v * v;
      const fast = HG.probeHit(th, k, L, objs);
      // brute force: march with step 1e-3, stop at the first sample inside any ball
      let bid = -1, bell = L;
      const ds = 1e-3;
      outer: for (let s = 0; s <= L; s += ds) {
        const p = HG.probe(th, k, s);
        for (let i = 0; i < objs.length; i++) {
          const o = objs[i];
          if (HG.gauge(HG.mul(HG.inv([o.x, o.y, o.t]), p)) < o.r) { bid = i; bell = s; break outer; }
        }
      }
      if (bid >= 0) hitsBrute++;
      if (bid !== fast.id) {
        // distinguish true misses from grazes: how deep does the brute-force winner (or the fast one) go?
        const depthOf = (id) => {
          if (id < 0) return 0;
          const o = objs[id]; let dmin = Infinity;
          for (let s = 0; s <= L; s += 5e-3) dmin = Math.min(dmin, HG.gauge(HG.mul(HG.inv([o.x, o.y, o.t]), HG.probe(th, k, s))) / o.r);
          return 1 - dmin;
        };
        const depth = Math.max(depthOf(bid), depthOf(fast.id));
        if (depth < 0.02) grazeMismatch++; else { idMismatch++; if (bad.length < 3) bad.push({ th, k, bid, bell, fid: fast.id, fell: fast.ell, depth }); }
      } else if (bid >= 0) {
        maxEll = Math.max(maxEll, Math.abs(bell - fast.ell));
      }
    }
  }
  console.log(`      eye hit test: ${probes} probes over 60 random scenes, ${hitsBrute} brute-force hits, ${idMismatch} real mismatches, ${grazeMismatch} grazes (penetration < 2% of radius) disagreeing`);
  if (bad.length) console.log('      mismatch examples:', JSON.stringify(bad));
  report('eye hit length vs brute-force march (step 1e-3), on agreeing hits', maxEll, 2e-3);
  report('eye hit test real mismatches (penetration >= 2% of radius)', idMismatch, 0);
}

// 8. inverse of the eye: preimages(probe(theta, kappa, l)) contains (theta, kappa, l), and every preimage it returns
//    really reaches the point. Counted against a brute-force scan of (u, l) for completeness.
{
  const L = 30, KMAX = 2;
  let recoverErr = 0, reachErr = 0, missing = 0, cases = 0, countMismatch = 0;
  seed = 4242;
  const wrapD = (a) => Math.abs(((a + 3 * Math.PI) % (2 * Math.PI)) - Math.PI);
  for (let n = 0; n < 3000; n++) {
    const th = U(-Math.PI, Math.PI), k = U(-KMAX, KMAX) * (rnd() < 0.2 ? 0.01 : 1), ell = U(0.05, L);
    const q = HG.probe(th, k, ell);
    if (Math.hypot(q[0], q[1]) < 1e-3) continue;          // on the t axis the preimage is a whole circle
    cases++;
    const pre = HG.preimages(q, L, KMAX);
    let bestD = Infinity;
    for (const p of pre) {
      bestD = Math.min(bestD, Math.max(wrapD(p.theta - th), Math.abs(p.kappa - k), Math.abs(p.ell - ell)));
      reachErr = Math.max(reachErr, dist(HG.probe(p.theta, p.kappa, p.ell), q) / (1 + Math.abs(q[2])));
    }
    if (bestD > 1e-6) { missing++; if (missing < 4) console.log('      not recovered:', JSON.stringify({ th, k, ell, pre })); }
    else recoverErr = Math.max(recoverErr, bestD);
    // brute force count: sign changes of t(u) - q_t along the ray, over u with the right r
    if (n < 300) {
      const r = Math.hypot(q[0], q[1]);
      let cnt = 0, prevG = null;
      for (let i = 1; i <= 400000; i++) {
        const u = (i / 400000) * KMAX * L * Math.sign(q[2] || 1) ;
        const s = 2 * Math.sin(u / 2) / u, l = r / Math.abs(s);
        const g = (l <= L && Math.abs(u / l) <= KMAX) ? HG.probe(0, u / l, l)[2] - q[2] : null;
        if (g !== null && prevG !== null && (g === 0 || (g > 0) !== (prevG > 0))) cnt++;
        prevG = g;
      }
      // the u < 0 side can only add the minimizing branch when q_t has the other sign, which the ray above already covers
      if (Math.abs(cnt - pre.length) > 0 && Math.abs(q[2]) > 1e-6) { countMismatch++; if (countMismatch < 4) console.log('      count:', cnt, 'vs', pre.length, JSON.stringify(q)); }
    }
  }
  report(`preimages recover the probe that made the point (${cases} points)`, recoverErr + missing, 1e-6);
  report('every preimage reaches the point, |probe - q| / (1 + |t|)', reachErr, 1e-8);
  report('number of preimages matches a brute-force root count along u (300 points)', countMismatch, 0);
}

console.log(failures ? `\n${failures} FAILED` : '\nall passed');
process.exit(failures ? 1 : 0);
