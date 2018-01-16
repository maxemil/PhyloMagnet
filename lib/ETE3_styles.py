from ete3 import *

def layout_node_color(node):
    """
    This layout set the color of a node if it has a "color" feature
    :param node:
    :return:
    """
    if node.is_leaf():
        name = node.name
        if "color" in node.features:
            N = faces.TextFace(name, fgcolor=node.color, fsize=12)
        else:
            N = faces.TextFace(name, fgcolor="#000000", fsize=12)
        faces.add_face_to_node(N, node, 0)


def tree_style_basic(layout=None, title=None):
    """
    It returns a tree style. It can take a layout function as argument
    :param layout: layout function
    :return: ETE3 TreeStyle Object
    """

    # Create an empty TreeStyle
    ts = TreeStyle()

    # Set title if provided
    if title is not None:
        ts.title.add_face(TextFace(title, fsize=20), column=2)

    # Set our custom layout function
    if layout is not None:
        ts.layout_fn = layout

    # Attributes
    # ts.mode = "c" # circular tree
    ts.scale = 400  # 400 pixels per branch length unit
    ts.branch_vertical_margin = 5  # 5 pixels between adjacent branches
    ts.show_leaf_name = False
    # ts.show_branch_length = True
    ts.show_branch_support = True
    ts.complete_branch_lines_when_necessary = False
    # ts.optimal_scale_level = "full"
    return ts


def node_style_basic():
    """
    It generates a NodeStyle
    :return: ETE3 NodeStyle Object
    """
    style = NodeStyle()
    style["fgcolor"] = "#000000"
    style["size"] = 0
    style["vt_line_color"] = "#000000"
    style["hz_line_color"] = "#000000"
    style["vt_line_width"] = 1
    style["hz_line_width"] = 1
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    return style
