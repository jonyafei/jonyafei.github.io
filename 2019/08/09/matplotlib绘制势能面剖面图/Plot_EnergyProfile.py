# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import xlrd
import sys
import re
from adjustText import adjust_text


def curve_points(y_small, y_large, direction):
    """产生余弦曲线插点"""
    if direction == "up":
        # 上升余弦曲线，产生纵轴y的插值点
        x_up = np.linspace(np.pi, 2 * np.pi, 50)
        y_curve = [y_small + (y_large - y_small) * (j + 1) / 2 for j in np.cos(x_up).tolist()]
    else:
        # 下降余弦曲线，产生纵轴y的插值点
        x_down = np.linspace(0, np.pi, 50)  # 下降余弦曲线横坐标
        y_curve = [y_small + (y_large - y_small) * (j + 1) / 2 for j in np.cos(x_down).tolist()]
    return y_curve


def interpolate_cos(x, y):
    """在反应路径驻点上产生一系列余弦曲线插值点"""
    x_new = []
    y_smooth = []

    # 1. 延伸起始点
    x_pre_temp = np.linspace(x[0] - 1, x[0], 50).tolist()
    # 下降余弦曲线，产生纵轴y的插值点，与初始两个点间曲线对称
    y_pre_temp = curve_points(y[0], y[1], "down")
    x_new = x_new + x_pre_temp
    y_smooth = y_smooth + y_pre_temp

    # 2. 中间点插值
    for i in range(len(x) - 1):
        x_new_temp = np.linspace(x[i], x[i + 1], 50).tolist()  # 产生横轴x的插值点
        if y[i] < y[i + 1]:
            y_smooth_temp = curve_points(y[i], y[i+1], "up")  # 上升余弦曲线，产生纵轴y的插值点
        else:
            y_smooth_temp = curve_points(y[i+1], y[i], "down")  # 下降余弦曲线，产生纵轴y的插值点
        x_new = x_new + x_new_temp  # 包含所有横坐标点的列表
        y_smooth = y_smooth + y_smooth_temp  # 包含所有纵坐标点的列表

    # 3. 延伸末端点
    x_post_temp = np.linspace(x[-1], x[-1] + 1, 50).tolist()
    # 上升余弦曲线，产生纵轴y的插值点，与最后两个点间曲线对称
    y_post_temp = curve_points(y[-1], y[-2], "up")
    x_new = x_new + x_post_temp
    y_smooth = y_smooth + y_post_temp

    return x_new, y_smooth  # 返回所有插值点的坐标


def plot_curve(x, y, color, path_label, TextLabel, FontSize=22):
    """绘制多条平滑曲线型能垒图"""
    x_new_array = []
    y_smooth_array = []
    for i in range(len(y)):  # 遍历所有列的能量值
        x_strip = []
        y_strip = []
        for j in range(len(y[i])):
            if y[i][j] != "":  # 剔除空数据
                y_strip.append(y[i][j])
                x_strip.append(x[j])

        plt.scatter(x_strip, y_strip, linewidth=6, color=color[i])  # 绘制散点图
        x_new, y_smooth = interpolate_cos(x_strip, y_strip)  # 调用插点函数，生成插点坐标
        plt.plot(x_new, y_smooth, linewidth=4, label=path_label[i], color=color[i])  # 绘制平滑曲线
        x_new_array = np.append(x_new_array, x_new)
        y_smooth_array = np.append(y_smooth_array, y_smooth)

        # 标记能量值，偏移量视具体情况而定
        if TextLabel == 'True':
            for j in range(len(x_strip)):
                texts = [plt.text(x_strip[j], y_strip[j], "{:.1f}".format(y_strip[j]), fontsize=FontSize,
                                  color=color[i])]  # 标记能量值，偏移量视具体情况而定
                adjust_text(texts, x_new_array, y_smooth_array)  # 调用adjust_text，尽可能避免文字文字、文字与曲线重叠


def line_split(x, y, color, TextLabel, FontSize=22):
    """绘制单条分段实线图"""
    y_new = []
    x_new = []
    # 1.生成新的XY坐标点，个数加倍
    for i in range(len(y)):
        if y[i] != "":  # 剔除空数据
            y_new.append(y[i])
            y_new.append(y[i])
            x_new.append(2*i+1)
            x_new.append(2*i+2)
    # 2.绘制实线折线图
    i = 0
    while i < len(y_new):
        x_line = [x_new[i], x_new[i+1]]
        y_line = [y_new[i], y_new[i+1]]
        plt.plot(x_line, y_line, linestyle='-', linewidth=6, color=color)
        i += 2
    # 3.添加能量值文本标签
    if TextLabel == 'True':
        for j in range(len(x)):
            if y[j] != "":
                plt.text(x[j] * 2 - 0.9, y[j] + 0.4, "{:.1f}".format(y[j]), fontsize=FontSize, color=color)
    return x_new, y_new


def plot_line_dot(x, y, color, path_label, TextLabel, FontSize=22):
    """绘制多条虚实折线图"""
    y_max, y_min = y_extreme(y)    # 获取y值的最大值和最小值
    y_bias = (y_max - y_min) / 50  # 获取文本标签y方向偏移量
    if isinstance(path_label, list):  # 多条路径的情况
        for i in range(len(y)):  # 遍历所有列的能量值
            # 绘制分段实线折线图
            x_new, y_new = line_split(x[i], y[i], color[i], TextLabel="False")
            # 绘制虚线折线图
            plt.plot(x_new, y_new, linestyle='--', linewidth=5, color=color[i], label=path_label[i])
            # 标记能量值，偏移量视具体情况而定
            if TextLabel == 'True':
                for j in range(len(x)):
                    if y[i][j] != "":
                        plt.text(x[j] * 2 - 0.9, y[i][j] + y_bias, "{:.1f}".format(y[i][j]), fontsize=FontSize,
                                 color=color[i])
    else:  # 单条路径的情况
        # 绘制分段实线折线图
        x_new, y_new = line_split(x, y, color, TextLabel="False")
        # 绘制虚线折线图
        plt.plot(x_new, y_new, linestyle='--', linewidth=5, color=color, label=path_label)

        # 标记能量值，偏移量视具体情况而定
        if TextLabel == 'True':
            for j in range(len(x)):
                if y[j] != "":
                    plt.text(x[j] * 2 - 0.9, y[j] + y_bias, "{:.1f}".format(y[j]), fontsize=FontSize, color=color)


def plot_line_curve(y_ini, _xtick_labels, color, PathLabel, TextLabel, FontSize=22):
    """单条曲线：根据中间体及过渡态类型绘制能垒图，对中间体绘制横线，对过渡态绘制曲线"""
    # 1. 数据预处理
    x_ini = [i * 2 + 2 for i in range(len(y_ini))]  # 产生x轴坐标
    x = []
    y = []
    for i in range(len(y_ini)):
        if y_ini[i] != "":
            if re.match(r"^TS", _xtick_labels[i]) is not None:  # 判断是否为TS数据点，若是，将y值添加到新的列表中
                y.append(y_ini[i])
                x.append(x_ini[i])
            else:  # 若否，将y值分两次添加到新的列表中
                y.append(y_ini[i])
                y.append(y_ini[i])
                x.append(x_ini[i]-0.5)
                x.append(x_ini[i]+0.5)

    # 2. 产生插值点
    x_new = []
    y_smooth = []
    for i in range(len(x) - 1):
        # 产生横轴x的插值点
        x_new_temp = np.linspace(x[i], x[i + 1], 50).tolist()
        if y[i] < y[i + 1]:
            # 上升余弦曲线，产生纵轴y的插值点
            y_smooth_temp = curve_points(y[i], y[i+1], "up")
        elif y[i] > y[i + 1]:
            # 下降余弦曲线，产生纵轴y的插值点
            y_smooth_temp = curve_points(y[i+1], y[i], "down")
        else:
            # 长横线
            y_smooth_temp = np.linspace(y[i], y[i+1], 50).tolist()
        x_new = x_new + x_new_temp  # 包含所有横坐标点的列表
        y_smooth = y_smooth + y_smooth_temp  # 包含所有纵坐标点的列表

    # 3. 绘制曲线
    plt.plot(x_new, y_smooth, linewidth=6, color=color, label=PathLabel)  # 绘制曲线

    # 4. 添加能量值文本标签
    if TextLabel == 'True':
        # 添加能量值文本标签
        for i in range(len(y_ini)):
            # 标记能量值，偏移量视具体情况而定
            if y_ini[i] != "":
                texts = [plt.text(x_ini[i] - 0.2, y_ini[i] + 0.06, "{:.1f}".format(y_ini[i]), fontsize=FontSize, color=color)]
        adjust_text(texts, x_new, y_smooth)  # 调用adjust_text，尽可能避免文字文字、文字与曲线重叠

    return x_ini  # 返回x标签点，用于绘制x轴标签


def plot_scatter(x_sticks, y, color, path_label, TextLabel, FontSize=22):
    """作散点图，并以长横线显示数据点"""
    for i in range(len(y)):  # 遍历所有列的能量值
        y_strip = []
        x_strip = []
        for j in range(len(y[i])):
            if y[i][j] != "":  # 剔除空数据
                x_strip.append(x_sticks[j])
                y_strip.append(y[i][j])
                # 添加文本标签
                if TextLabel == "True":
                    plt.text(x_sticks[j] - 0.3, y[i][j] + 0.06, "{:.1f}".format(y[i][j]), fontsize=FontSize, color=color[i])
        # 绘制其他自旋态，画横线
        plt.scatter(x_strip, y_strip, linewidth=6, color=color[i], label=path_label[i], marker='_', s=1200)


def y_extreme(y):
    """返回y列表中的最大值和最小值"""
    nest = "False"
    for i in y:
        if isinstance(i, list):
            nest = "True"  # 若是嵌套列表，给nest赋值为True
            break
    if nest == "True":
        temp = sum(y, [])  # 展开y列表
    else:
        temp = y
    temp = [i for i in temp if i != ""]  # 剔除空数据
    y_max, y_min = max(temp), min(temp)
    return y_max, y_min


def y_list_min(y):
    """寻找每行的最小值，并返回一个最小值的列表"""
    y_T = np.array(y).T.tolist()  # 转置y列表数据
    y_min = []
    for i in range(len(y_T)):
        for j in range(y_T[i].count("")):
            y_T[i].remove("")  # 删除空数据
        y_min.append(min(list(map(float, y_T[i]))))  # 返回每行的最小值，添加到y_min_list列表中
    return y_min


# 1. 导入数据
ExcelFile = xlrd.open_workbook(sys.argv[1])  # 读取Excel数据
sheet = ExcelFile.sheet_by_index(0)  # 读取Excel的第一个sheet
_xtick_labels = sheet.col_values(0)[5:]  # 读取第一列数据，反应路径驻点的名称
pic_title = sheet.row_values(0)[1]  # 读取图片标题
X_title = sheet.row_values(1)[1]  # 读取X轴标题
Y_title = sheet.row_values(2)[1]  # 读取Y轴标题
path_label = sheet.row_values(3)[1:]  # 读取第一行数据，不同反应路径的名称
color = sheet.row_values(4)[1:]
x = [i+1 for i in range(len(_xtick_labels))]  # 生成横坐标
y = []
for i in range(len(sheet.row_values(0))-1):
    y.append(sheet.col_values(i+1)[5:])  # 读取除第一列外的所有列数据，即纵坐标能量值

y_min_list = y_list_min(y)  # 返回所有中间体及过渡态中最稳定自旋态的能量值

# 2. 绘制图像
plt.figure(figsize=(15, 9), dpi=80)  # 设置图片大小及分辨率

plot_style = int(sys.argv[2])
TextLabel = str(sys.argv[3])  # 是否添加坐标点对应的数值文本标签
FontSize = 22  # 设置能量值文本大小，默认值为22，可根据需要修改
Axis_FontSize = 20  # 设置XY轴标签及标题大小，默认值为20，可根据需要修改
if plot_style == 1:  # 绘制平滑曲线
    plot_curve(x, y, color, path_label, TextLabel, FontSize)
    plt.xlim(x[0] - 0.5, x[-1] + 0.5)  # x轴刻度范围
    plt.xticks(x, _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
elif plot_style == 11:
    x_new, y_smooth = interpolate_cos([i * 2 - 0.5 for i in x], y_min_list)  # 调用插点函数，生成插点坐标
    plt.plot(x_new, y_smooth, color="grey", label=None, linewidth=6)  # 绘制平滑曲线
    plot_scatter([i * 2 - 0.5 for i in x], y, color, path_label, TextLabel, FontSize)
    plt.xlim(x[0] * 2 - 1.5, x[-1] * 2 + 1)  # x轴刻度范围
    plt.xticks([i * 2 - 0.5 for i in x], _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
elif plot_style == 2:  # 绘制虚实折线
    plot_line_dot(x, y, color, path_label, TextLabel, FontSize)
    plt.xlim(x[0] * 2 - 1.5, x[-1] * 2 + 1)  # x轴刻度范围
    plt.xticks([i * 2 - 0.5 for i in x], _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
elif plot_style == 22:  # # 对不同自旋态，只给最稳定态绘制虚实折线
    plot_line_dot(x, y_min_list, "grey", None, TextLabel="False")  #
    for i in range(len(y)):
        line_split(x, y[i], color[i], TextLabel="False")
    plot_scatter([i * 2 - 0.5 for i in x], y, color, path_label, TextLabel, FontSize)
    plt.xlim(x[0] * 2 - 1.5, x[-1] * 2 + 1)  # x轴刻度范围
    plt.xticks([i * 2 - 0.5 for i in x], _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
elif plot_style == 3:  # 绘制横线&平滑曲线
    for i in range(len(y)):  # 遍历所有列的能量值
        x_sticks = plot_line_curve(y[i], _xtick_labels, color[i], path_label[i], TextLabel, FontSize)
        plt.xticks(x_sticks, _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
else:  # 33 对不同自旋态，只能最稳定态绘制横线&平滑线，其他点绘制横线
    x_sticks = plot_line_curve(y_min_list, _xtick_labels, "grey", None, TextLabel="False")  # 绘制平滑线
    plt.xticks(x_sticks, _xtick_labels, fontsize=Axis_FontSize)  # x轴标签
    plot_scatter(x_sticks, y, color, path_label, TextLabel, FontSize)

# 若x轴标签过长产生重叠，可设置旋转角度，比如rotation=-90, HorizontalAlignment="right"

# 3. 图片设置
y_max, y_min = y_extreme(y)  # 获取y值的最大值和最小值
y_scale = (y_max - y_min) / 10  # y轴延伸长度
plt.ylim(y_min - y_scale, y_max + y_scale)  # y轴刻度范围
plt.yticks(fontsize=Axis_FontSize)  # y轴标签
# plt.xlabel(X_title, fontsize=Axis_FontSize)  # 横轴标题
plt.ylabel(Y_title, fontsize=Axis_FontSize)  # 纵轴标题
# plt.title(pic_title, fontsize=Axis_FontSize)  # 图标题

plt.legend(fontsize=Axis_FontSize-2, loc="upper right")  # 添加图例,位置在左上角
plt.tight_layout()  # 图像外部边缘的调整

# plt.show()  # 展示图片
plt.savefig("./EnergyProfile.png")  # 保存图片到当前目录
