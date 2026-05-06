#!/usr/bin/env python3
"""
Chemical Element Query Tool - 化学元素查询工具

查询元素周期表中化学元素的详细信息。
支持通过元素名称（中英文）、符号或原子序数查询。

Usage:
    python element_query.py --query Au
    python element_query.py --query 金
    python element_query.py --query Gold
    python element_query.py --query 79 --format markdown
"""

import argparse
import json
import sys
from typing import Union, Optional

# 元素周期表数据（1-118号元素）
ELEMENTS = {
    1: {"symbol": "H", "name": "Hydrogen", "name_cn": "氢", "mass": 1.008, "description": "氢是宇宙中最丰富的元素，约占宇宙质量的75%。它是恒星核聚变的主要燃料，也是生命体系中水（H₂O）和有机分子的基本组成部分。", "eleconfig": "1s1"},
    2: {"symbol": "He", "name": "Helium", "name_cn": "氦", "mass": 4.0026, "description": "氦是第二轻的元素，无色无味不可燃的惰性气体。常用于填充气球、飞艇，以及作为超导磁体的冷却剂（液氦）。", "eleconfig": "1s2"},
    3: {"symbol": "Li", "name": "Lithium", "name_cn": "锂", "mass": 6.94, "description": "锂是最轻的金属元素，质软，可用刀切割。广泛用于锂离子电池、精神疾病药物（碳酸锂）以及某些合金中。", "eleconfig": "[He] 2s1"},
    4: {"symbol": "Be", "name": "Beryllium", "name_cn": "铍", "mass": 9.0122, "description": "铍是一种轻质、坚硬的金属，具有优异的刚度和导热性。用于航空航天结构材料、X射线窗口和核反应堆中子减速剂。", "eleconfig": "[He] 2s2"},
    5: {"symbol": "B", "name": "Boron", "name_cn": "硼", "mass": 10.81, "description": "硼是一种类金属元素，自然界中以硼酸盐形式存在。用于硼硅玻璃（耐热玻璃）、洗涤剂、半导体掺杂和某些植物营养素。", "eleconfig": "[He] 2s2 2p1"},
    6: {"symbol": "C", "name": "Carbon", "name_cn": "碳", "mass": 12.011, "description": "碳是生命的基础元素，能形成复杂多样的有机分子。存在同素异形体：金刚石、石墨、石墨烯、富勒烯和碳纳米管。", "eleconfig": "[He] 2s2 2p2"},
    7: {"symbol": "N", "name": "Nitrogen", "name_cn": "氮", "mass": 14.007, "description": "氮占地球大气的78%，是所有蛋白质和核酸的关键组成元素。液氮广泛用于低温冷却和食品冷冻保存。", "eleconfig": "[He] 2s2 2p3"},
    8: {"symbol": "O", "name": "Oxygen", "name_cn": "氧", "mass": 15.999, "description": "氧是地壳中含量最丰富的元素，也是水和大多数生物体的核心组成。约占大气的21%，是呼吸和燃烧过程的必需元素。", "eleconfig": "[He] 2s2 2p4"},
    9: {"symbol": "F", "name": "Fluorine", "name_cn": "氟", "mass": 18.998, "description": "氟是最活泼的非金属元素，能与几乎所有元素反应。用于制造聚四氟乙烯（Teflon）、牙膏添加剂（氟化物）和制冷剂。", "eleconfig": "[He] 2s2 2p5"},
    10: {"symbol": "Ne", "name": "Neon", "name_cn": "氖", "mass": 20.180, "description": "氖是一种无色单原子分子惰性气体，通电时发出红橙色光，广泛用于霓虹灯和广告灯管。", "eleconfig": "[He] 2s2 2p6"},
    11: {"symbol": "Na", "name": "Sodium", "name_cn": "钠", "mass": 22.990, "description": "钠是碱金属元素，质软活泼，与水剧烈反应。氯化钠（食盐）是常见的化合物。钠离子在神经传导和体液平衡中起关键作用。", "eleconfig": "[Ne] 3s1"},
    12: {"symbol": "Mg", "name": "Magnesium", "name_cn": "镁", "mass": 24.305, "description": "镁是轻质金属，燃烧时发出耀眼的白光。叶绿素的核心金属离子就是镁。广泛用于合金、阻燃剂和医药（泻盐）。", "eleconfig": "[Ne] 3s2"},
    13: {"symbol": "Al", "name": "Aluminium", "name_cn": "铝", "mass": 26.982, "description": "铝是地壳中含量最丰富的金属，质轻耐腐蚀。广泛用于航空航天、建筑材料、包装（铝箔和易拉罐）和输电线。", "eleconfig": "[Ne] 3s2 3p1"},
    14: {"symbol": "Si", "name": "Silicon", "name_cn": "硅", "mass": 28.085, "description": "硅是地壳中含量第二丰富的元素，半导体工业的基石。硅芯片是现代电子设备和计算机的核心。也用于玻璃和陶瓷制造。", "eleconfig": "[Ne] 3s2 3p2"},
    15: {"symbol": "P", "name": "Phosphorus", "name_cn": "磷", "mass": 30.974, "description": "磷是生命必需元素，存在于DNA、RNA、ATP和骨骼中。存在白磷、红磷、黑磷等同素异形体。广泛用于肥料和洗涤剂。", "eleconfig": "[Ne] 3s2 3p3"},
    16: {"symbol": "S", "name": "Sulfur", "name_cn": "硫", "mass": 32.06, "description": "硫是黄色固体非金属，有特殊臭味。存在于蛋白质和维生素中。用于制造硫酸（最重要的工业化学品之一）、橡胶硫化等。", "eleconfig": "[Ne] 3s2 3p4"},
    17: {"symbol": "Cl", "name": "Chlorine", "name_cn": "氯", "mass": 35.45, "description": "氯是黄绿色有毒气体，强氧化剂。用于水消毒、漂白剂、聚氯乙烯（PVC）塑料和多种化学品制造。", "eleconfig": "[Ne] 3s2 3p5"},
    18: {"symbol": "Ar", "name": "Argon", "name_cn": "氩", "mass": 39.948, "description": "氩是空气中含量最多的惰性气体（约0.93%）。用于填充白炽灯泡、焊接保护和半导体制造中的惰性气氛。", "eleconfig": "[Ne] 3s2 3p6"},
    19: {"symbol": "K", "name": "Potassium", "name_cn": "钾", "mass": 39.098, "description": "钾是碱金属，对生物体至关重要。钾离子维持细胞渗透压和神经肌肉兴奋性。广泛用于肥料和某些药物。", "eleconfig": "[Ar] 4s1"},
    20: {"symbol": "Ca", "name": "Calcium", "name_cn": "钙", "mass": 40.078, "description": "钙是第五丰富的地壳元素，骨骼和牙齿的主要成分（羟基磷灰石）。也用于水泥、石灰和金属还原剂。", "eleconfig": "[Ar] 4s2"},
    21: {"symbol": "Sc", "name": "Scandium", "name_cn": "钪", "mass": 44.956, "description": "钪是过渡金属，质轻强度高。用于航空航天合金、高强度放电灯和某些激光材料。", "eleconfig": "[Ar] 3d1 4s2"},
    22: {"symbol": "Ti", "name": "Titanium", "name_cn": "钛", "mass": 47.867, "description": "钛是强度高、耐腐蚀、生物相容性好的金属。广泛用于航空航天部件、人工关节、牙科植入物和化工设备。", "eleconfig": "[Ar] 3d2 4s2"},
    23: {"symbol": "V", "name": "Vanadium", "name_cn": "钒", "mass": 50.942, "description": "钒是高强度钢的合金元素。五氧化二钒是硫酸生产的重要催化剂。也用于储能电池（钒液流电池）。", "eleconfig": "[Ar] 3d3 4s2"},
    24: {"symbol": "Cr", "name": "Chromium", "name_cn": "铬", "mass": 51.996, "description": "铬是坚硬耐腐蚀的金属，不锈钢的关键合金成分。镀铬提供高光泽和耐磨表面。六价铬有毒，三价铬是人体必需微量元素。", "eleconfig": "[Ar] 3d5 4s1"},
    25: {"symbol": "Mn", "name": "Manganese", "name_cn": "锰", "mass": 54.938, "description": "锰是重要的工业金属，用于钢铁生产（脱氧脱硫）和干电池（二氧化锰）。也是人体必需的微量元素。", "eleconfig": "[Ar] 3d5 4s2"},
    26: {"symbol": "Fe", "name": "Iron", "name_cn": "铁", "mass": 55.845, "description": "铁是使用最广泛的金属，人类文明的基础。血红蛋白的核心，负责氧气运输。不锈钢和各类合金的主要成分。", "eleconfig": "[Ar] 3d6 4s2"},
    27: {"symbol": "Co", "name": "Cobalt", "name_cn": "钴", "mass": 58.933, "description": "钴是磁性金属，用于高性能合金、永磁体和锂离子电池正极材料（钴酸锂）。维生素B12的核心元素。", "eleconfig": "[Ar] 3d7 4s2"},
    28: {"symbol": "Ni", "name": "Nickel", "name_cn": "镍", "mass": 58.693, "description": "镍是不锈钢的主要合金元素，也用于电池（镍氢、镍镉、镍氢）、电镀和催化剂。", "eleconfig": "[Ar] 3d8 4s2"},
    29: {"symbol": "Cu", "name": "Copper", "name_cn": "铜", "mass": 63.546, "description": "铜是人类最早使用的金属之一，优良的导电导热材料。广泛用于电线、管道、硬币和合金（青铜、黄铜）。", "eleconfig": "[Ar] 3d10 4s1"},
    30: {"symbol": "Zn", "name": "Zinc", "name_cn": "锌", "mass": 65.38, "description": "锌是活泼金属，用于镀锌钢板防锈、干电池负极和合金（黄铜）。是人体必需的微量元素，参与免疫功能。", "eleconfig": "[Ar] 3d10 4s2"},
    31: {"symbol": "Ga", "name": "Gallium", "name_cn": "镓", "mass": 69.723, "description": "镓是熔点很低（29.76°C）的金属，在手中即可熔化。用于半导体（砷化镓）、LED和太阳能电池。", "eleconfig": "[Ar] 3d10 4s2 4p1"},
    32: {"symbol": "Ge", "name": "Germanium", "name_cn": "锗", "mass": 72.63, "description": "锗是重要的半导体材料，也用于光纤通信系统和红外光学器件。锗晶体管是早期电子工业的基础。", "eleconfig": "[Ar] 3d10 4s2 4p2"},
    33: {"symbol": "As", "name": "Arsenic", "name_cn": "砷", "mass": 74.922, "description": "砷是著名的有毒元素，但微量砷对某些生物是必需的。用于半导体（砷化镓）、木材防腐剂和某些杀虫剂。", "eleconfig": "[Ar] 3d10 4s2 4p3"},
    34: {"symbol": "Se", "name": "Selenium", "name_cn": "硒", "mass": 78.96, "description": "硒是人体必需的微量元素，抗氧化酶谷胱甘肽过氧化物酶的组成成分。用于光电池、复印机和玻璃着色。", "eleconfig": "[Ar] 3d10 4s2 4p4"},
    35: {"symbol": "Br", "name": "Bromine", "name_cn": "溴", "mass": 79.904, "description": "溴是唯一在常温下呈液态的非金属元素，红棕色液体。用于阻燃剂、消毒剂、摄影化学品和药物制造。", "eleconfig": "[Ar] 3d10 4s2 4p5"},
    36: {"symbol": "Kr", "name": "Krypton", "name_cn": "氪", "mass": 83.798, "description": "氪是惰性气体，通电时发出明亮的白色光。用于高性能灯泡、摄影闪光灯和某些激光器。", "eleconfig": "[Ar] 3d10 4s2 4p6"},
    37: {"symbol": "Rb", "name": "Rubidium", "name_cn": "铷", "mass": 85.468, "description": "铷是碱金属，非常活泼。用于原子钟、光电池和某些特种玻璃。", "eleconfig": "[Kr] 5s1"},
    38: {"symbol": "Sr", "name": "Strontium", "name_cn": "锶", "mass": 87.62, "description": "锶是碱土金属，燃烧时发出红色火焰。用于焰火、信号弹、磁铁和某些陶瓷材料。", "eleconfig": "[Kr] 5s2"},
    39: {"symbol": "Y", "name": "Yttrium", "name_cn": "钇", "mass": 88.906, "description": "钇是稀土金属，用于LED荧光粉、激光晶体（YAG）和高温超导材料。", "eleconfig": "[Kr] 4d1 5s2"},
    40: {"symbol": "Zr", "name": "Zirconium", "name_cn": "锆", "mass": 91.224, "description": "锆是耐腐蚀金属，用于核反应堆包壳材料、耐火材料和陶瓷（氧化锆）。", "eleconfig": "[Kr] 4d2 5s2"},
    41: {"symbol": "Nb", "name": "Niobium", "name_cn": "铌", "mass": 92.906, "description": "铌是超导金属，用于超导磁体（如MRI）、航空航天合金和电子器件。", "eleconfig": "[Kr] 4d4 5s1"},
    42: {"symbol": "Mo", "name": "Molybdenum", "name_cn": "钼", "mass": 95.95, "description": "钼是高熔点金属，人体必需的微量元素（钼酶）。用于高强度钢和催化剂。", "eleconfig": "[Kr] 4d5 5s1"},
    43: {"symbol": "Tc", "name": "Technetium", "name_cn": "锝", "mass": 98, "description": "锝是人造元素，自然界中几乎不存在。用于医学影像（⁹⁹ᵐTc放射性示踪剂），是核医学中使用最广泛的同位素。", "eleconfig": "[Kr] 4d5 5s2"},
    44: {"symbol": "Ru", "name": "Ruthenium", "name_cn": "钌", "mass": 101.07, "description": "钌是铂族金属，用于耐磨电触点、太阳能电池和某些催化剂。", "eleconfig": "[Kr] 4d7 5s1"},
    45: {"symbol": "Rh", "name": "Rhodium", "name_cn": "铑", "mass": 102.91, "description": "铑是铂族贵金属，用于汽车催化转化器、珠宝电镀和化学催化剂。", "eleconfig": "[Kr] 4d8 5s1"},
    46: {"symbol": "Pd", "name": "Palladium", "name_cn": "钯", "mass": 106.42, "description": "钯是铂族金属，能吸收大量氢气。用于汽车催化转化器、牙科材料和电子元件。", "eleconfig": "[Kr] 4d10"},
    47: {"symbol": "Ag", "name": "Silver", "name_cn": "银", "mass": 107.87, "description": "银是导电导热性最好的金属，抗菌材料。用于珠宝、货币、电子触点和摄影胶片。", "eleconfig": "[Kr] 4d10 5s1"},
    48: {"symbol": "Cd", "name": "Cadmium", "name_cn": "镉", "mass": 112.41, "description": "镉是有毒重金属。用于镍镉电池、颜料（镉黄/镉红）和某些合金。", "eleconfig": "[Kr] 4d10 5s2"},
    49: {"symbol": "In", "name": "Indium", "name_cn": "铟", "mass": 114.82, "description": "铟是柔软的银白色金属。用于制造ITO透明导电膜（触摸屏、液晶显示器）和某些半导体。", "eleconfig": "[Kr] 4d10 5s2 5p1"},
    50: {"symbol": "Sn", "name": "Tin", "name_cn": "锡", "mass": 118.71, "description": "锡是人类最早使用的金属之一。用于焊接、马口铁（镀锡钢板）和青铜合金。", "eleconfig": "[Kr] 4d10 5s2 5p2"},
    51: {"symbol": "Sb", "name": "Antimony", "name_cn": "锑", "mass": 121.76, "description": "锑是脆性金属，用于铅酸电池合金、阻燃剂和半导体材料（锑化铟）。", "eleconfig": "[Kr] 4d10 5s2 5p3"},
    52: {"symbol": "Te", "name": "Tellurium", "name_cn": "碲", "mass": 127.60, "description": "碲是准金属，用于太阳能电池（碲化镉）、热电材料和某些合金。", "eleconfig": "[Kr] 4d10 5s2 5p4"},
    53: {"symbol": "I", "name": "Iodine", "name_cn": "碘", "mass": 126.90, "description": "碘是紫黑色固体，升华产生紫色蒸气。人体必需微量元素（甲状腺激素）。用于消毒剂（碘酒）和摄影化学品。", "eleconfig": "[Kr] 4d10 5s2 5p5"},
    54: {"symbol": "Xe", "name": "Xenon", "name_cn": "氙", "mass": 131.29, "description": "氙是最重的稳定惰性气体。用于高强度放电灯（氙气灯）、麻醉剂和离子推进器。", "eleconfig": "[Kr] 4d10 5s2 5p6"},
    55: {"symbol": "Cs", "name": "Cesium", "name_cn": "铯", "mass": 132.91, "description": "铯是最活泼的碱金属，熔点最低（28.5°C）。用于原子钟、光电管和石油钻探液。", "eleconfig": "[Xe] 6s1"},
    56: {"symbol": "Ba", "name": "Barium", "name_cn": "钡", "mass": 137.33, "description": "钡是碱土金属，燃烧发绿色光。用于烟火、信号弹、钡餐（X射线造影剂）和石油钻探。", "eleconfig": "[Xe] 6s2"},
    57: {"symbol": "La", "name": "Lanthanum", "name_cn": "镧", "mass": 138.91, "description": "镧是镧系元素的第一个，用于镍氢电池、光学玻璃和催化剂（石油精炼）。", "eleconfig": "[Xe] 5d1 6s2"},
    58: {"symbol": "Ce", "name": "Cerium", "name_cn": "铈", "mass": 140.12, "description": "铈是最丰富的稀土元素。用于打火石（铈铁合金）、抛光粉、催化剂和玻璃着色。", "eleconfig": "[Xe] 4f1 5d1 6s2"},
    59: {"symbol": "Pr", "name": "Praseodymium", "name_cn": "镨", "mass": 140.91, "description": "镨是稀土金属，用于制造高强度永磁体、航空合金和黄色玻璃颜料。", "eleconfig": "[Xe] 4f3 6s2"},
    60: {"symbol": "Nd", "name": "Neodymium", "name_cn": "钕", "mass": 144.24, "description": "钕是制造最强永磁体（钕铁硼）的关键元素。也用于激光晶体和玻璃着色。", "eleconfig": "[Xe] 4f4 6s2"},
    61: {"symbol": "Pm", "name": "Promethium", "name_cn": "钷", "mass": 145, "description": "钷是放射性稀土元素，自然界中不存在。用于核电池、荧光材料和某些研究应用。", "eleconfig": "[Xe] 4f5 6s2"},
    62: {"symbol": "Sm", "name": "Samarium", "name_cn": "钐", "mass": 150.36, "description": "钐是稀土金属，用于永磁体、核反应堆控制棒和某些催化剂。", "eleconfig": "[Xe] 4f6 6s2"},
    63: {"symbol": "Eu", "name": "Europium", "name_cn": "铕", "mass": 151.96, "description": "铕是稀土元素，用于红色荧光粉（电视和节能灯）、核反应堆控制材料。", "eleconfig": "[Xe] 4f7 6s2"},
    64: {"symbol": "Gd", "name": "Gadolinium", "name_cn": "钆", "mass": 157.25, "description": "钆是稀土金属，用于MRI造影剂、核反应堆控制棒和磁制冷材料。", "eleconfig": "[Xe] 4f7 5d1 6s2"},
    65: {"symbol": "Tb", "name": "Terbium", "name_cn": "铽", "mass": 158.93, "description": "铽是稀土金属，用于荧光粉、磁致伸缩材料（声纳）和某些激光器。", "eleconfig": "[Xe] 4f9 6s2"},
    66: {"symbol": "Dy", "name": "Dysprosium", "name_cn": "镝", "mass": 162.50, "description": "镝是稀土金属，用于高性能永磁体、核反应堆控制棒和磁致伸缩材料。", "eleconfig": "[Xe] 4f10 6s2"},
    67: {"symbol": "Ho", "name": "Holmium", "name_cn": "钬", "mass": 164.93, "description": "钬是稀土金属，具有最高的磁矩。用于激光器、核反应堆控制棒和某些磁体。", "eleconfig": "[Xe] 4f11 6s2"},
    68: {"symbol": "Er", "name": "Erbium", "name_cn": "铒", "mass": 167.26, "description": "铒是稀土金属，用于光纤放大器（EDFA，长距离光通信）、激光器和核应用。", "eleconfig": "[Xe] 4f12 6s2"},
    69: {"symbol": "Tm", "name": "Thulium", "name_cn": "铥", "mass": 168.93, "description": "铥是稀土金属，用于便携式X射线设备、激光器和某些核应用。", "eleconfig": "[Xe] 4f13 6s2"},
    70: {"symbol": "Yb", "name": "Ytterbium", "name_cn": "镱", "mass": 173.05, "description": "镱是稀土金属，用于光纤放大器、原子钟和某些激光器。", "eleconfig": "[Xe] 4f14 6s2"},
    71: {"symbol": "Lu", "name": "Lutetium", "name_cn": "镥", "mass": 174.97, "description": "镥是最重的镧系元素，用于石油精炼催化剂和某些核应用。", "eleconfig": "[Xe] 4f14 5d1 6s2"},
    72: {"symbol": "Hf", "name": "Hafnium", "name_cn": "铪", "mass": 178.49, "description": "铪是与锆非常相似的金属，用于核反应堆控制棒和某些特种合金。", "eleconfig": "[Xe] 4f14 5d2 6s2"},
    73: {"symbol": "Ta", "name": "Tantalum", "name_cn": "钽", "mass": 180.95, "description": "钽是耐腐蚀金属，用于外科手术植入物、电容器和高温合金。", "eleconfig": "[Xe] 4f14 5d3 6s2"},
    74: {"symbol": "W", "name": "Tungsten", "name_cn": "钨", "mass": 183.84, "description": "钨是熔点最高的金属（3422°C）。用于白炽灯丝、硬质合金刀具和高炉耐火材料。", "eleconfig": "[Xe] 4f14 5d4 6s2"},
    75: {"symbol": "Re", "name": "Rhenium", "name_cn": "铼", "mass": 186.21, "description": "铼是稀有过渡金属，用于喷气发动机高温合金和石油精炼催化剂。", "eleconfig": "[Xe] 4f14 5d5 6s2"},
    76: {"symbol": "Os", "name": "Osmium", "name_cn": "锇", "mass": 190.23, "description": "锇是密度最大的天然元素（22.59 g/cm³）。用于合金硬化剂和某些电触点。", "eleconfig": "[Xe] 4f14 5d6 6s2"},
    77: {"symbol": "Ir", "name": "Iridium", "name_cn": "铱", "mass": 192.22, "description": "铱是最耐腐蚀的金属之一。用于火花塞、高温 crucibles 和某些电触点。", "eleconfig": "[Xe] 4f14 5d7 6s2"},
    78: {"symbol": "Pt", "name": "Platinum", "name_cn": "铂", "mass": 195.08, "description": "铂是贵金属，化学性质极稳定。用于汽车催化转化器、珠宝、电极和某些抗癌药物（顺铂）。", "eleconfig": "[Xe] 4f14 5d9 6s1"},
    79: {"symbol": "Au", "name": "Gold", "name_cn": "金", "mass": 196.97, "description": "金是贵金属，自古以来用作货币和珠宝。优良的导电性和抗腐蚀性，用于电子工业、牙科和航天器涂层。", "eleconfig": "[Xe] 4f14 5d10 6s1"},
    80: {"symbol": "Hg", "name": "Mercury", "name_cn": "汞", "mass": 200.59, "description": "汞是唯一在常温下呈液态的金属。用于温度计、荧光灯和某些化学过程。有毒，需谨慎处理。", "eleconfig": "[Xe] 4f14 5d10 6s2"},
    81: {"symbol": "Tl", "name": "Thallium", "name_cn": "铊", "mass": 204.38, "description": "铊是有毒重金属。曾用于鼠药和杀虫剂，现主要用于电子和光学玻璃。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p1"},
    82: {"symbol": "Pb", "name": "Lead", "name_cn": "铅", "mass": 207.2, "description": "铅是重金属，密度大、柔软。用于蓄电池、辐射屏蔽和某些合金。有毒，需避免摄入。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p2"},
    83: {"symbol": "Bi", "name": "Bismuth", "name_cn": "铋", "mass": 208.98, "description": "铋是重金属中毒性最低的。用于胃药（碱式水杨酸铋/Pepto-Bismol）、化妆品和低熔点合金。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p3"},
    84: {"symbol": "Po", "name": "Polonium", "name_cn": "钋", "mass": 209, "description": "钋是放射性元素，由居里夫妇发现。用于太空探测器热源（核电池）和某些工业测量设备。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p4"},
    85: {"symbol": "At", "name": "Astatine", "name_cn": "砹", "mass": 210, "description": "砹是放射性卤素，自然界中极稀有。用于某些癌症放射性治疗研究。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p5"},
    86: {"symbol": "Rn", "name": "Radon", "name_cn": "氡", "mass": 222, "description": "氡是放射性惰性气体，铀衰变产物。是室内空气污染源之一，与肺癌风险相关。", "eleconfig": "[Xe] 4f14 5d10 6s2 6p6"},
    87: {"symbol": "Fr", "name": "Francium", "name_cn": "钫", "mass": 223, "description": "钫是最活泼的碱金属，极不稳定，半衰期仅22分钟。用于某些原子物理研究。", "eleconfig": "[Rn] 7s1"},
    88: {"symbol": "Ra", "name": "Radium", "name_cn": "镭", "mass": 226, "description": "镭是放射性碱土金属，居里夫妇发现。曾用于发光涂料（表盘），现主要用于癌症放疗。", "eleconfig": "[Rn] 7s2"},
    89: {"symbol": "Ac", "name": "Actinium", "name_cn": "锕", "mass": 227, "description": "锕是锕系元素的第一个，强放射性。用于某些癌症放射性治疗研究。", "eleconfig": "[Rn] 6d1 7s2"},
    90: {"symbol": "Th", "name": "Thorium", "name_cn": "钍", "mass": 232.04, "description": "钍是放射性金属，潜在核燃料（钍基反应堆）。也用于合金和光学玻璃。", "eleconfig": "[Rn] 6d2 7s2"},
    91: {"symbol": "Pa", "name": "Protactinium", "name_cn": "镤", "mass": 231.04, "description": "镤是放射性金属，极稀有。主要用于科学研究。", "eleconfig": "[Rn] 5f2 6d1 7s2"},
    92: {"symbol": "U", "name": "Uranium", "name_cn": "铀", "mass": 238.03, "description": "铀是主要核燃料，用于核电站和核武器。自然界中含量比银和汞还多。", "eleconfig": "[Rn] 5f3 6d1 7s2"},
    93: {"symbol": "Np", "name": "Neptunium", "name_cn": "镎", "mass": 237, "description": "镎是第一个人造超铀元素。用于某些核研究和探测器。", "eleconfig": "[Rn] 5f4 6d1 7s2"},
    94: {"symbol": "Pu", "name": "Plutonium", "name_cn": "钚", "mass": 244, "description": "钚是重要核材料，用于核武器和核反应堆（如旅行者号探测器核电池）。", "eleconfig": "[Rn] 5f6 7s2"},
    95: {"symbol": "Am", "name": "Americium", "name_cn": "镅", "mass": 243, "description": "镅是人造元素，用于烟雾探测器（镅-241）和某些密度测量设备。", "eleconfig": "[Rn] 5f7 7s2"},
    96: {"symbol": "Cm", "name": "Curium", "name_cn": "锔", "mass": 247, "description": "锔是人造元素，为纪念居里夫妇命名。用于某些太空探测器热源和科学研究。", "eleconfig": "[Rn] 5f7 6d1 7s2"},
    97: {"symbol": "Bk", "name": "Berkelium", "name_cn": "锫", "mass": 247, "description": "锫是人造元素，用于某些基础科学研究。", "eleconfig": "[Rn] 5f9 7s2"},
    98: {"symbol": "Cf", "name": "Californium", "name_cn": "锎", "mass": 251, "description": "锎是人造元素，某些同位素强中子源。用于中子活化分析、石油勘探和癌症放疗研究。", "eleconfig": "[Rn] 5f10 7s2"},
    99: {"symbol": "Es", "name": "Einsteinium", "name_cn": "锿", "mass": 252, "description": "锿是人造元素，为纪念爱因斯坦命名。仅用于科学研究。", "eleconfig": "[Rn] 5f11 7s2"},
    100: {"symbol": "Fm", "name": "Fermium", "name_cn": "镄", "mass": 257, "description": "镄是人造元素，为纪念费米命名。仅用于科学研究。", "eleconfig": "[Rn] 5f12 7s2"},
    101: {"symbol": "Md", "name": "Mendelevium", "name_cn": "钔", "mass": 258, "description": "钔是人造元素，为纪念门捷列夫命名。仅用于科学研究。", "eleconfig": "[Rn] 5f13 7s2"},
    102: {"symbol": "No", "name": "Nobelium", "name_cn": "锘", "mass": 259, "description": "锘是人造元素，为纪念诺贝尔命名。仅用于科学研究。", "eleconfig": "[Rn] 5f14 7s2"},
    103: {"symbol": "Lr", "name": "Lawrencium", "name_cn": "铹", "mass": 262, "description": "铹是人造元素，为纪念劳伦斯命名。仅用于科学研究。", "eleconfig": "[Rn] 5f14 7s2 7p1"},
    104: {"symbol": "Rf", "name": "Rutherfordium", "name_cn": "𬬻", "mass": 267, "description": "𬬻是超重元素，为纪念卢瑟福命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d2 7s2"},
    105: {"symbol": "Db", "name": "Dubnium", "name_cn": "𬭊", "mass": 268, "description": "𬭊是超重元素，为纪念杜布纳联合核研究所命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d3 7s2"},
    106: {"symbol": "Sg", "name": "Seaborgium", "name_cn": "𬭳", "mass": 269, "description": "𬭳是超重元素，为纪念西博格命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d4 7s2"},
    107: {"symbol": "Bh", "name": "Bohrium", "name_cn": "𬭛", "mass": 270, "description": "𬭛是超重元素，为纪念玻尔命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d5 7s2"},
    108: {"symbol": "Hs", "name": "Hassium", "name_cn": "𬭶", "mass": 269, "description": "𬭶是超重元素，为纪念黑森州命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d6 7s2"},
    109: {"symbol": "Mt", "name": "Meitnerium", "name_cn": "鿏", "mass": 278, "description": "鿏是超重元素，为纪念迈特纳命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d7 7s2"},
    110: {"symbol": "Ds", "name": "Darmstadtium", "name_cn": "𫟼", "mass": 281, "description": "𫟼是超重元素，为纪念达姆施塔特命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d8 7s2"},
    111: {"symbol": "Rg", "name": "Roentgenium", "name_cn": "𬬭", "mass": 282, "description": "𬬭是超重元素，为纪念伦琴命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d9 7s2"},
    112: {"symbol": "Cn", "name": "Copernicium", "name_cn": "鿔", "mass": 285, "description": "鿔是超重元素，为纪念哥白尼命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2"},
    113: {"symbol": "Nh", "name": "Nihonium", "name_cn": "鿭", "mass": 286, "description": "鿭是超重元素，由日本RIKEN研究所发现，名称意为'日本'。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p1"},
    114: {"symbol": "Fl", "name": "Flerovium", "name_cn": "𫓧", "mass": 289, "description": "𫓧是超重元素，为纪念弗廖罗夫命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p2"},
    115: {"symbol": "Mc", "name": "Moscovium", "name_cn": "镆", "mass": 289, "description": "镆是超重元素，为纪念莫斯科州命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p3"},
    116: {"symbol": "Lv", "name": "Livermorium", "name_cn": "𫟷", "mass": 293, "description": "𫟷是超重元素，为纪念劳伦斯利弗莫尔国家实验室命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p4"},
    117: {"symbol": "Ts", "name": "Tennessine", "name_cn": "鿬", "mass": 294, "description": "鿬是超重元素，为纪念田纳西州命名。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p5"},
    118: {"symbol": "Og", "name": "Oganesson", "name_cn": "鿫", "mass": 294, "description": "鿫是最重的已知元素，为纪念奥加涅相命名。理论预测可能不是惰性气体。仅用于核物理研究。", "eleconfig": "[Rn] 5f14 6d10 7s2 7p6"},
}


# 构建查找索引
_SYMBOL_TO_NUMBER = {data["symbol"]: num for num, data in ELEMENTS.items()}
_NAME_TO_NUMBER = {data["name"].lower(): num for num, data in ELEMENTS.items()}
_NAME_CN_TO_NUMBER = {data["name_cn"]: num for num, data in ELEMENTS.items()}


def _expand_eleconfig(config: str) -> str:
    """将简写电子排布展开为完整形式"""
    noble_gases = {
        "[He]": "1s2",
        "[Ne]": "1s2 2s2 2p6",
        "[Ar]": "1s2 2s2 2p6 3s2 3p6",
        "[Kr]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
        "[Xe]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6",
        "[Rn]": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6 6s2 4f14 5d10 6p6",
    }
    
    for noble, expansion in noble_gases.items():
        if config.startswith(noble):
            return config.replace(noble, expansion, 1)
    return config


def _parse_eleconfig_to_dict(config: str) -> dict:
    """将电子排布字符串解析为结构化字典"""
    expanded = _expand_eleconfig(config)
    eleconfig_dict = {}
    
    # 解析格式如 "1s2 2s2 2p6"
    parts = expanded.split()
    for part in parts:
        # 提取轨道和电子数，如 "3d10" -> 轨道="3d", 电子数=10
        import re
        match = re.match(r'(\d+[spdf])(\d+)', part)
        if match:
            orbital = match.group(1)
            electrons = int(match.group(2))
            eleconfig_dict[orbital] = electrons
    
    return eleconfig_dict


def query_element(query: Union[str, int]) -> dict:
    """
    查询化学元素信息。
    
    Args:
        query: 元素查询条件，可以是：
            - 元素符号（如 "Au", "H"）
            - 中文名称（如 "金", "氢"）
            - 英文名称（如 "Gold", "Hydrogen"）
            - 原子序数（如 79, 1）
    
    Returns:
        dict: 元素详细信息，包含原子序数、符号、名称、原子量、描述、电子排布等
    """
    atomic_number = None
    
    # 尝试作为整数解析（原子序数）
    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        atomic_number = int(query)
    elif isinstance(query, str):
        query_str = query.strip()
        
        # 尝试作为符号查找
        if query_str in _SYMBOL_TO_NUMBER:
            atomic_number = _SYMBOL_TO_NUMBER[query_str]
        # 尝试作为中文名称查找
        elif query_str in _NAME_CN_TO_NUMBER:
            atomic_number = _NAME_CN_TO_NUMBER[query_str]
        # 尝试作为英文名称查找（不区分大小写）
        elif query_str.lower() in _NAME_TO_NUMBER:
            atomic_number = _NAME_TO_NUMBER[query_str.lower()]
    
    # 返回结果
    if atomic_number and atomic_number in ELEMENTS:
        data = ELEMENTS[atomic_number].copy()
        data["number"] = atomic_number
        
        # 展开电子排布为完整形式
        data["eleconfig"] = _expand_eleconfig(data["eleconfig"])
        
        # 生成结构化电子排布字典
        data["eleconfig_dict"] = _parse_eleconfig_to_dict(data["eleconfig"])
        
        return {
            "Status": "success",
            "Message": "Element found",
            **data
        }
    else:
        return {
            "Status": "error",
            "Message": f"Element not found: {query}. Supported queries: symbol (e.g., Au), Chinese name (e.g., 金), English name (e.g., Gold), or atomic number (e.g., 79)."
        }


def format_markdown(result: dict) -> str:
    """将结果格式化为Markdown表格"""
    if result["Status"] == "error":
        return f"❌ {result['Message']}"
    
    lines = [
        f"**Element Information: {result['name']} ({result['name_cn']})**",
        "",
        "| Property | Value |",
        "|----------|-------|",
        f"| **Atomic Number** | {result['number']} |",
        f"| **Symbol** | {result['symbol']} |",
        f"| **Name** | {result['name']} ({result['name_cn']}) |",
        f"| **Atomic Mass** | {result['mass']} amu |",
        f"| **Electron Configuration** | `{result['eleconfig']}` |",
        "",
        "**Description:**",
        result['description'],
        "",
        "**Electron Configuration (Structured):**",
        "```json",
        json.dumps(result['eleconfig_dict'], indent=2, ensure_ascii=False),
        "```",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Chemical Element Query Tool - Query periodic table element information'
    )
    parser.add_argument('--query', '-q', required=True,
                        help='Element query: symbol (Au), name (Gold/金), or atomic number (79)')
    parser.add_argument('--format', '-f', choices=['json', 'markdown'], default='json',
                        help='Output format (default: json)')
    
    args = parser.parse_args()
    
    # 查询元素
    result = query_element(args.query)
    
    # 格式化输出
    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)
    
    print(output)
    
    return 0 if result["Status"] == "success" else 1


if __name__ == '__main__':
    sys.exit(main())
