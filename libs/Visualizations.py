from sage.all import *
from math import sin, cos, pi, sqrt


# ============================================================
# 基础工具函数
# ============================================================

def edge_to_uvl(e):
    """
    兼容 Sage 中边可能是 (u,v) 或 (u,v,label) 的情况。
    """
    if len(e) == 2:
        return e[0], e[1], None
    else:
        return e[0], e[1], e[2]


def normalize_edge_key(u, v, lab=None):
    """
    无向边的标准 key，用于查找 edge_graphs。
    """
    a, b = sorted([u, v], key=str)
    return (a, b, lab)


def shift_segment(p, q, shift):
    """
    把线段 pq 沿法向平移 shift。
    用于画多重边时错开边。
    """
    x1, y1 = p
    x2, y2 = q

    dx = x2 - x1
    dy = y2 - y1
    length = sqrt(dx * dx + dy * dy)

    if length == 0:
        return p, q

    nx = -dy / length
    ny = dx / length

    return (
        (x1 + shift * nx, y1 + shift * ny),
        (x2 + shift * nx, y2 + shift * ny)
    )




def bezier_control_point(p, q, curvature):
    """
    给端点 p, q 构造一个二次 Bezier 曲线的控制点。
    curvature > 0 / < 0 决定弯曲方向和程度。
    """
    x1, y1 = p
    x2, y2 = q

    mx = (x1 + x2) / 2
    my = (y1 + y2) / 2

    dx = x2 - x1
    dy = y2 - y1
    d = sqrt(dx*dx + dy*dy)

    if d == 0:
        return (mx, my)

    # 单位法向量
    nx = -dy / d
    ny = dx / d

    # 控制点沿法向偏移
    cx = mx + curvature * d * nx
    cy = my + curvature * d * ny

    return (cx, cy)


def quadratic_bezier_point(p, c, q, t=0.5):
    """
    二次 Bezier 曲线在参数 t 处的点。
    """
    x1, y1 = p
    cx, cy = c
    x2, y2 = q

    x = (1-t)**2 * x1 + 2*(1-t)*t * cx + t**2 * x2
    y = (1-t)**2 * y1 + 2*(1-t)*t * cy + t**2 * y2

    return (x, y)


def draw_bezier_edge(p, q, curvature=0.0, color="blue", thickness=1.2, zorder=10):
    """
    画一条二次 Bezier 边。
    curvature=0 时退化成直线。
    返回: Graphics, control_point
    """
    P = Graphics()

    if abs(curvature) < 1e-12:
        P += line([p, q], color=color, thickness=thickness, zorder=zorder)
        c = ((p[0]+q[0])/2, (p[1]+q[1])/2)
        return P, c

    c = bezier_control_point(p, q, curvature)

    # Sage 里 bezier_path 的输入是一个 path 列表，
    # 每段 path 用 [起点, 控制点, 终点]
    P += bezier_path([[p, c, q]], color=color, thickness=thickness, zorder=zorder)

    return P, c


def multiedge_curvatures(m, base=0.22):
    """
    给 m 条平行边分配一组对称的弯曲参数。

    例子:
    m=1 -> [0]
    m=2 -> [-0.11, 0.11]
    m=3 -> [-0.22, 0, 0.22]
    m=4 -> [-0.33, -0.11, 0.11, 0.33]
    """
    if m == 1:
        return [0.0]

    return [base * (k - (m - 1) / 2) for k in range(m)]




def edge_to_uvl(e):
    if len(e) == 2:
        return e[0], e[1], None
    return e[0], e[1], e[2]


def shift_segment(p, q, shift):
    x1, y1 = p
    x2, y2 = q

    dx = x2 - x1
    dy = y2 - y1
    d = sqrt(dx*dx + dy*dy)

    if d == 0:
        return p, q

    nx = -dy / d
    ny = dx / d

    return (
        (x1 + shift * nx, y1 + shift * ny),
        (x2 + shift * nx, y2 + shift * ny)
    )


def trim_segment(p, q, r):
    """
    把边的两端截掉一点，避免被顶点白圆完全盖住。
    """
    x1, y1 = p
    x2, y2 = q

    dx = x2 - x1
    dy = y2 - y1
    d = sqrt(dx*dx + dy*dy)

    if d <= 2*r:
        return p, q

    ux = dx / d
    uy = dy / d

    return (
        (x1 + r * ux, y1 + r * uy),
        (x2 - r * ux, y2 - r * uy)
    )


def local_circle_pos(
    H,
    center,
    radius,
    inner_scale=0.55,
    vertex_layout="circle",
    padding=0.20
):
    """
    小图顶点布局。

    vertex_layout:
        "circle"  : 保留原来的手动圆形布局，默认值；
        "sage"    : 使用 H.get_pos() 或 H.layout()；
        "spring"  : 使用 Sage 的 spring layout；
        "planar"  : 使用 Sage 的 planar layout；
        "circular": 使用 Sage 的 circular layout；
        dict      : 直接使用你传入的位置字典。
    """
    V = list(H.vertices(sort=False))
    n = len(V)

    cx, cy = center

    if n == 0:
        return {}

    if n == 1:
        return {V[0]: (cx, cy)}

    # --------------------------------------------------
    # 情况 1：保留原来的圆形布局
    # --------------------------------------------------
    if vertex_layout == "circle":
        r = inner_scale * radius

        return {
            v: (
                cx + r * cos(pi / 2 + 2*pi*i/n),
                cy + r * sin(pi / 2 + 2*pi*i/n)
            )
            for i, v in enumerate(V)
        }

    # --------------------------------------------------
    # 情况 2：用户直接传位置字典
    # --------------------------------------------------
    if isinstance(vertex_layout, dict):
        raw_pos = vertex_layout

    # --------------------------------------------------
    # 情况 3：使用 Sage 的布局
    # --------------------------------------------------
    else:
        raw_pos = None

        # "sage" 时，优先尊重 H 已有的位置
        if vertex_layout == "sage":
            try:
                raw_pos = H.get_pos()
            except Exception:
                raw_pos = None

            if raw_pos is None or len(raw_pos) == 0:
                try:
                    raw_pos = H.layout()
                except Exception:
                    raw_pos = H.layout_circular()

        # "spring", "planar", "circular" 等
        else:
            try:
                raw_pos = H.layout(layout=vertex_layout)
            except Exception:
                try:
                    raw_pos = H.layout(vertex_layout)
                except Exception:
                    try:
                        raw_pos = H.layout()
                    except Exception:
                        raw_pos = H.layout_circular()

    # --------------------------------------------------
    # 把 Sage 原始布局缩放到小圆内部
    # --------------------------------------------------
    raw_pos = {v: raw_pos[v] for v in V}

    xs = [float(raw_pos[v][0]) for v in V]
    ys = [float(raw_pos[v][1]) for v in V]

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2

    width = xmax - xmin
    height = ymax - ymin

    denom = max(width, height, 1e-9)

    usable_radius = (1 - padding) * inner_scale * radius

    # bbox 宽度要放入直径，所以乘 2
    scale = 2 * usable_radius / denom

    return {
        v: (
            cx + scale * (float(raw_pos[v][0]) - xmid),
            cy + scale * (float(raw_pos[v][1]) - ymid)
        )
        for v in V
    }


def draw_small_graph(
    H,
    center,
    radius,
    node_radius=None,
    edge_thickness=0.45,
    vertex_outline_thickness=0.35,
    vertex_labels=False,
    label_fontsize=7,
    inner_scale=0.55,
    vertex_layout="circle",
    debug_edges=False,
):
    """
    手动画小图。

    画图顺序：
    1. 小图内部边；
    2. 小图顶点白色圆盘；
    3. 小图顶点黑色边框；
    4. 可选顶点标签。

    所有点和边都使用 local_circle_pos 给出的同一套坐标。
    """
    P = Graphics()

    if H is None:
        return P

    V = list(H.vertices(sort=False))
    E = list(H.edges(sort=False))

    if len(V) == 0:
        return P

    pos = local_circle_pos(
        H,
        center=center,
        radius=radius,
        inner_scale=inner_scale,
        vertex_layout=vertex_layout
    )

    if node_radius is None:
        node_radius = 0.055 * radius

    # 按无向边分组，处理小图内部多重边
    edge_groups = {}

    for e in E:
        u, v, lab = edge_to_uvl(e)
        key = tuple(sorted([u, v], key=str))
        edge_groups.setdefault(key, []).append((u, v, lab))

    # 先画小图内部边
    for key, edges in edge_groups.items():
        m = len(edges)

        for k, (u, v, lab) in enumerate(edges):
            p = pos[u]
            q = pos[v]

            if u == v:
                x, y = p
                # loop_r = 2.2 * node_radius
                loop_r = 1.4 * node_radius

                P += circle(
                    (x + loop_r, y + loop_r),
                    loop_r,
                    fill=False,
                    color="red" if debug_edges else "black",
                    thickness=edge_thickness,
                    zorder=50
                )

            else:
                curvs = multiedge_curvatures(m, base=0.35)

                for (u, v, lab), curv in zip(edges, curvs):
                    p = pos[u]
                    q = pos[v]

                    if u == v:
                        x, y = p
                        loop_r = 2.2 * node_radius
                        P += circle(
                            (x + loop_r, y + loop_r),
                            loop_r,
                            fill=False,
                            color="red" if debug_edges else "black",
                            thickness=edge_thickness,
                            zorder=50
                        )
                    else:
                        # 先把边的起终点缩进一点，避免被顶点盖住
                        # p_trim, q_trim = trim_segment(p, q, 1.15 * node_radius)
                        p_trim, q_trim = p, q

                        G_edge, ctrl = draw_bezier_edge(
                            p_trim, q_trim,
                            curvature=curv,
                            color="red" if debug_edges else "black",
                            thickness=edge_thickness,
                            zorder=50
                        )
                        P += G_edge

    # 再画小图顶点
    for v in V:
        x, y = pos[v]

        P += circle(
            (x, y),
            0.5* node_radius,
            fill=True,
            color="white",
            zorder=60
        )

        P += circle(
            (x, y),
            node_radius,
            fill=False,
            color="black",
            thickness=vertex_outline_thickness,
            zorder=61
        )

    if vertex_labels:
        for v in V:
            x, y = pos[v]

            P += text(
                str(v),
                (x, y),
                fontsize=label_fontsize,
                color="black",
                zorder=70
            )

    return P

def edge_length(p, q):
    return sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)


def rescale_outer_pos_by_min_edge_length(Big, pos, target_edge_length):
    """
    把外层顶点坐标整体放大，使最短外层边长度至少为 target_edge_length。
    """
    lengths = []

    seen = set()

    for e in Big.edges(sort=False):
        u, v, lab = edge_to_uvl(e)

        if u == v:
            continue

        key = tuple(sorted([u, v], key=str))
        if key in seen:
            continue
        seen.add(key)

        if u in pos and v in pos:
            d = edge_length(pos[u], pos[v])
            if d > 1e-9:
                lengths.append(d)

    if len(lengths) == 0:
        return pos

    current_min = min(lengths)

    scale = target_edge_length / current_min

    xs = [float(p[0]) for p in pos.values()]
    ys = [float(p[1]) for p in pos.values()]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)

    return {
        v: (
            cx + scale * (float(p[0]) - cx),
            cy + scale * (float(p[1]) - cy)
        )
        for v, p in pos.items()
    }

def auto_figsize_from_pos(
    pos,
    vertex_radius=0.5,
    edge_graph_radius=0.25,
    units_per_inch=0.9,
    min_width=8,
    min_height=5,
    margin=1.0
):
    """
    根据外层图坐标范围自动生成 figsize。
    
    units_per_inch 越小，图越大、越清楚。
    推荐 0.7 ~ 1.0。
    """
    xs = [float(p[0]) for p in pos.values()]
    ys = [float(p[1]) for p in pos.values()]

    r = max(vertex_radius, edge_graph_radius)

    width_units = max(xs) - min(xs) + 2 * (r + margin)
    height_units = max(ys) - min(ys) + 2 * (r + margin)

    fig_w = max(min_width, width_units / units_per_inch)
    fig_h = max(min_height, height_units / units_per_inch)

    return [fig_w, fig_h]


# ============================================================
# 画嵌套多重图
# ============================================================

def rescale_outer_pos(pos, scale=1.0):
    """
    把外层大图顶点位置相对于中心整体放大。
    scale > 1 会增长外层边长。
    """
    if scale == 1.0:
        return pos

    xs = [float(p[0]) for p in pos.values()]
    ys = [float(p[1]) for p in pos.values()]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)

    return {
        v: (
            cx + scale * (float(p[0]) - cx),
            cy + scale * (float(p[1]) - cy)
        )
        for v, p in pos.items()
    }

def draw_nested_multigraph(
    Big,
    vertex_graphs,
    edge_graphs=None,
    big_pos=None,
    vertex_radius=0.35,
    edge_graph_radius=0.16,
    vertex_node_radius=None,
    edge_node_radius=None,
    show_outer_vertex_labels=True,
    outer_label_fontsize=10,
    debug=False,
    figsize=8,
    outer_scale=1.0,
    target_edge_length=None,
    outer_edge_thickness=0.7,
    outer_vertex_outline_thickness=1.0,
    edge_vertex_outline_thickness=0.8,
    small_vertex_edge_thickness=0.45,
    small_edge_edge_thickness=0.35,
    small_vertex_outline_thickness=0.35,
    small_edge_outline_thickness=0.30,
):
    """
    Big:
        外层 Sage multigraph。

    vertex_graphs:
        dict，形如
        {
            外层顶点: 该顶点里面要画的小图
        }

    edge_graphs:
        dict，形如
        {
            (u, v, label): 该边上要画的小图
        }

        注意：由于 Big 是无向图，所以 (u,v,label) 和 (v,u,label) 都可以。
        程序内部会自动标准化。

    big_pos:
        外层大图顶点位置。
        None 时自动调用 Big.layout()。
    """
    if edge_graphs is None:
        edge_graphs = {}

    # 规范化 edge_graphs 的 key
    normalized_edge_graphs = {}

    for key, H in edge_graphs.items():
        if len(key) == 2:
            u, v = key
            lab = None
        else:
            u, v, lab = key

        normalized_edge_graphs[normalize_edge_key(u, v, lab)] = H

    if big_pos is None:
        big_pos = Big.layout()

    big_pos = rescale_outer_pos(big_pos, scale=outer_scale)

    if target_edge_length is not None:
        big_pos = rescale_outer_pos_by_min_edge_length(
        Big,
        big_pos,
        target_edge_length=target_edge_length
    )

    if vertex_node_radius is None:
        vertex_node_radius = 0.055 * vertex_radius

    if edge_node_radius is None:
        edge_node_radius = 0.060 * edge_graph_radius

    P = Graphics()

    # --------------------------------------------------------
    # debug 检查：确认小图对象本身有没有边
    # --------------------------------------------------------
    if debug:
        print("===== vertex_graphs check =====")
        for v, H in vertex_graphs.items():
            print(
                "outer vertex =", v,
                "small graph vertices =", H.order(),
                "small graph edges =", H.size(),
                "edges =", H.edges(sort=False)
            )

        print("===== edge_graphs check =====")
        for k, H in normalized_edge_graphs.items():
            print(
                "outer edge =", k,
                "small graph vertices =", H.order(),
                "small graph edges =", H.size(),
                "edges =", H.edges(sort=False)
            )

    # --------------------------------------------------------
    # 外层边分组，用于处理多重边
    # --------------------------------------------------------
    edge_groups = {}

    for e in Big.edges(sort=False):
        u, v, lab = edge_to_uvl(e)
        key = normalize_edge_key(u, v, None)
        edge_groups.setdefault(key, []).append((u, v, lab))

    # --------------------------------------------------------
    # 先画外层大图的边
    # --------------------------------------------------------
    edge_draw_data = []

    for key, edges in edge_groups.items():
        m = len(edges)
        curvs = multiedge_curvatures(m, base=0.5) #0.22

        for (u, v, lab), curv in zip(edges, curvs):
            p = big_pos[u]
            q = big_pos[v]

            G_edge, ctrl = draw_bezier_edge(
                p, q,
                curvature=curv,
                color="blue",
                thickness=outer_edge_thickness,
                zorder=5
            )
            P += G_edge

            # 用 Bezier 曲线中点作为“边上小图”的位置
            mid = quadratic_bezier_point(p, ctrl, q, t=0.5)

            edge_draw_data.append((u, v, lab, p, q, ctrl, mid, curv))

    # --------------------------------------------------------
    # 再画边上的小图
    # --------------------------------------------------------
    for u, v, lab, p, q, ctrl, mid, curv in edge_draw_data:
        mx, my = mid

        edge_key = normalize_edge_key(u, v, lab)
        H_edge = normalized_edge_graphs.get(edge_key, None)

        if H_edge is not None:
            # 边上小图的白色背景圆
            P += circle(
                (mx, my),
                edge_graph_radius,
                fill=True,
                color="white"
            )

            P += circle(
                (mx, my),
                edge_graph_radius,
                fill=False,
                color="blue",
                thickness=edge_vertex_outline_thickness
            )

            P += draw_small_graph(
                H_edge,
                center=(mx, my),
                radius=edge_graph_radius,
                node_radius=0.035 * edge_graph_radius,
                edge_thickness=small_edge_edge_thickness,
                vertex_outline_thickness=small_edge_outline_thickness,
                inner_scale=0.62,
                vertex_layout="spring",
                vertex_labels=False,
            )

    # --------------------------------------------------------
    # 最后画外层顶点的大圆和里面的小图
    # --------------------------------------------------------
    for v in Big.vertices(sort=False):
        x, y = big_pos[v]

        # 外层顶点的大白圆，必须在小图之前画
        P += circle(
            (x, y),
            vertex_radius,
            fill=True,
            color="white"
        )

        P += circle(
            (x, y),
            vertex_radius,
            fill=False,
            color="blue",
            thickness=outer_vertex_outline_thickness
        )

        H_v = vertex_graphs.get(v, None)

        if H_v is not None:
            P += draw_small_graph(
                H_v,
                center=(x, y),
                radius=vertex_radius,
                node_radius=0.035 * vertex_radius,
                edge_thickness=small_vertex_edge_thickness,
                vertex_outline_thickness=small_vertex_outline_thickness,
                inner_scale=0.62,
                vertex_layout="spring",
                vertex_labels=False,
            )

        if show_outer_vertex_labels:
            P += text(
                str(v),
                (x, y - 1.25 * vertex_radius),
                color="blue",
                fontsize=outer_label_fontsize
            )

    P.axes(False)
    P.set_aspect_ratio(1)
    P.show(figsize=figsize)

    return P