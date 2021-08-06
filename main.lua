local err = 1.0E-10;

local function isNaN(n)
    return n ~= n;
end

local function solveLinear(a, b)
    return -b / a;
end

local function solveQuadratic(a, b, c)
    local k = -b / (2 * a);
    local u2 = k * k - c / a;
    if u2 > -err and u2 < err then
        return;
    else
        local u = u2 ^ 0.5;
        local r1, r2 = k - u, k + u;
        return r1, r2;
    end
end

local function solveCubic(a, b, c, d)
    local k = -b / (3 * a);
	local p = (3 * a * c - b * b) / (9 * a * a);
	local q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (54 * a * a * a);
	local r = p * p * p + q * q;
	local s = r ^ 0.5 + q;
    if s > -err and s < err then
        if q < 0 then
            return k + (-2 * q) ^ (1 / 3);
        else
            return k - (2 * q) ^ (1 / 3);
        end
    elseif r < 0 then
        local m = (-p) ^ 0.5
		local d = math.atan2((-r) ^ 0.5, q) / 3;
		local u = m * math.cos(d);
		local v = m * math.sin(d);
		return k - 2 * u, k + u - 1.7320508075688772 * v, k + u + 1.7320508075688772 * v;
    elseif s < 0 then
        local m = -(-s) ^ (1 / 3);
		return k + p / m - m;
    else
        local m = s ^ (1 / 3);
		return k + p / m - m;
    end
end

local function solveQuartic(a, b, c, d, e)
    local k = -b / (4 * a);
	local p = (8 * a * c - 3 * b * b) / (8 * a * a);
	local q = (b * b * b + 8 * a * a * d - 4 * a * b * c) / (8 * a * a * a);
	local r = (16 * a * a * b * b * c + 256 * a * a * a * a * e - 3 * a * b * b * b * b - 64 * a * a * a * b * d) / (256 * a * a * a * a * a);
    local h0, h1, h2 = solveCubic(1, 2 * p, p * p - 4 * r, -q * q);
    local s = h2 or h0;
    if s < err then
        local f0, f1 = solveQuadratic(1, p, r);
        if not f1 or f1 < 0 then
            return;
        else
            local f = f1 ^ 0.5;
            return k - f, k + f;
        end
    else
        local h = s ^ 0.5;
		local f = (h * h * h + h * p - q) / (2 * h);
		if f > -err and f < err then
			return k - h, k;
        else
            local r0, r1 = solveQuadratic(1, h, f);
            local r2, r3 = solveQuadratic(1, -h, r / f);
            if r0 and r2 then
                return k + r0, k + r1, k + r2, k + r3;
            elseif r0 then
                return k + r0, k + r1;
            elseif r2 then
                return k + r2, k + r3;
            else
                return;
            end
        end
    end
end

local function trajectory(origin, target, gravity, projectileSpeed)
    --origin : Vector3
    --target : Vector3
    --target : number(sign = positive)
    --projectileSpeed : number(sign = positive)
    local ox, oy, oz = origin.x, origin.y, origin.z;
    local tx, ty, tz = target.x, target.y, target.z;
    local tvx, tvy, tvz = tx - ox, ty - oy, tz - oz;
    local h = tx - ox;
    local j = tz - oz;
    local k = ty - oy;
    local l = -0.5 * gravity;

    local c0 = l * l;
    local c1 = -2 * tvy * l;
    local c2 = tvy*tvy - 2*k*l - projectileSpeed*projectileSpeed + tvx*tvx + tvz*tvz;
    local c3 = 2*k*tvy + 2*h*tvx + 2*j*tvz;
    local c4 = k*k + h*h + j*j;

    local t0, t1, t2, t3 = solveQuartic(c0, c1, c2, c3, c4);
    local d, e, f;
    if t0 and t0 >= 0 then
        d = ((h+tvx*t0)/t0);
        e = ((k+tvy*t0-l*t0*t0)/t0);
        f = ((j+tvz*t0)/t0);
        return Vector3.new(d, e, f), t0;
    end
    if t1 and t1 >= 0 then
        d = ((h+tvx*t1)/t1);
        e = ((k+tvy*t1-l*t1*t1)/t1);
        f = ((j+tvz*t1)/t1);
        return Vector3.new(d, e, f), t1;
    end
    if t2 and t2 >= 0 then
        d = ((h+tvx*t2)/t2);
        e = ((k+tvy*t2-l*t2*t2)/t2);
        f = ((j+tvz*t2)/t2);
        return Vector3.new(d, e, f), t2;
    end
    if t3 and t3 >= 0 then
        d = ((h+tvx*t3)/t3);
        e = ((k+tvy*t3-l*t3*t3)/t3);
        f = ((j+tvz*t3)/t3);
        return Vector3.new(d, e, f), t3;
    end
end


return {
    polynomials = {
        linear = solveLinear,
        quadratic = solveQuadratic,
        cubic = solveCubic,
        quartic = solveQuartic
    },
    projectile = {
        trajectory  = trajectory,

    }
};
