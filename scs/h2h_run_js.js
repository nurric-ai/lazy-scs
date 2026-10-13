// usage: node run_js.js <file>   (one string per line) -> prints length, lowerBound, seconds
const fs = require('fs');
const { solve } = require(process.argv[3] || './js/src/index.js');
const strs = fs.readFileSync(process.argv[2], 'utf8').split(/\r?\n/).filter(s => s.length);
const t = process.hrtime.bigint();
const r = solve(strs);
const secs = Number(process.hrtime.bigint() - t) / 1e9;
console.log(JSON.stringify({ length: r.superstring.length, lowerBound: r.lowerBound, fallback: r.fallback, seconds: secs, heapMB: process.memoryUsage().rss / 2**20 }));
