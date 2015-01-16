#!/usr/bin/env python

def print_tree(t):
  def p(node, d=1):
    if node is not t.nil:
      print "\t"*d + str(node.left)
      p(node.left, d+1)
      print "\t"*d + str(node.right)
      p(node.right, d+1)
  if t.root is not t.nil:
    print str(t.root)
    p(t.root)

class Color:
  def __repr__(self):
    return self.color

class Red(Color):
  def __init__(self):
    self.color = "red"

class Black(Color):
  def __init__(self):
    self.color = "black"

class RedBlack():
  red = Red()
  black = Black()
 
class RedBlackNode(RedBlack):
  def __init__(self, key=None, 
                     high=None,
                     max=None,
                     parent=None,
                     left=None,
                     right=None,
                     color=RedBlack.red):
    self.parent = parent
    self.left = left
    self.right = right
    self.key = key
    self.high = high
    self.max = max
    self.color = color
  
  def __repr__(self):
    return "(%s, %s)%s: %s" % (str(self.key), str(self.high), 
                               str(self.max), str(self.color))

class RedBlackTree():
  def __init__(self, nodetype=RedBlackNode):
    self.nil = nodetype(color=nodetype.black)
    self.root = self.nil
    self.root.parent = self.nil
    self.root.left = self.nil
    self.root.right = self.nil

  def left_rotate(self, x):
    y = x.right
    x.right = y.left
    if y.left is not self.nil:
      y.left.parent = x
    y.parent = x.parent
    if x.parent is self.nil:
      self.root = y
    elif x is x.parent.left:
      x.parent.left = y
    else:
      x.parent.right = y
    y.left = x
    x.parent = y
    x.max=max(x.high, x.left.max, x.right.max)
    y.max=max(y.high, y.left.max, y.right.max)

  def right_rotate(self, x):
    y = x.left
    x.left = y.right
    if y.right is not self.nil:
      y.right.parent = x
    y.parent = x.parent
    if x.parent is self.nil:
      self.root = y
    elif x is x.parent.right:
      x.parent.right = y
    else:
      x.parent.left = y
    y.right = x
    x.parent = y
    x.max=max(x.high, x.left.max, x.right.max)
    y.max=max(y.high, y.left.max, y.right.max)
  
  def insert(self, z):
    y = self.nil
    x = self.root
    while x is not self.nil:
      y = x
      if z.key < x.key:
        x = x.left
      else:
        x = x.right
    z.parent = y
    if y is self.nil:
      self.root = z
    elif z.key < y.key:
      y.left = z
    else:
      y.right = z
    z.left = self.nil
    z.right = self.nil
    z.color = RedBlackNode.red
    z.max = z.high
    p = z.parent
    while p is not self.nil:
      p.max = max(p.high, p.left.max, p.right.max)
      p = p.parent
    self.fixup(z)

  def fixup(self, z):
    red, black = RedBlackNode.red, RedBlackNode.black
    while z.parent.color == red:
      if z.parent is z.parent.parent.left:
        y = z.parent.parent.right
        if y.color == red:
          z.parent.color = black
          y.color = black
          z.parent.parent.color = red
          z = z.parent.parent
        else:
          if z is z.parent.right:
            z = z.parent
            self.left_rotate(z)
          z.parent.color = black
          z.parent.parent.color = red
          self.right_rotate(z.parent.parent)
      else:
        y = z.parent.parent.left
        if y.color == red:
          z.parent.color = black
          y.color = black
          z.parent.parent.color = red
          z = z.parent.parent
        else:
          if z is z.parent.left:
            z = z.parent
            self.right_rotate(z)
          z.parent.color = black
          z.parent.parent.color = red
          self.left_rotate(z.parent.parent)
    self.root.color = black
  
  def search(self, p):
    x = self.root
    while x is not self.nil and (p<x.key or p >=x.high):
      if x.left is not self.nil and x.left.max > p:
        x = x.left
      else:
        x = x.right
    return x
  
  def __contains__(self, k):
    if self.search(k) is not self.nil:
      return True
    else:
      return False

class IntervalSet():
  def __init__(self):
    self.t = RedBlackTree()
  
  def insert(self, low, high):
    self.t.insert(RedBlackNode(key=low, high=high))
  
  def search(self, x):
    n = self.t.search(x)
    if n is not self.t.nil:
      return (n.key, n.high)
    else:
      return None
  
  def __contains__(self, x):
    print_tree(self.t)
    return x in self.t

def read_vcf(f):
  res = {}
  for l in f:
    if l[0] != "#":
      spl= l.split()
      chr = spl[0]
      pos = int(spl[1])
      if res.has_key(chr):
        res[chr].append(pos-1)
      else:
        res[chr] = [pos-1]
  return res

def make_bed_intset(f):
  res = {}
  for l in f:
    spl = l.split()
    chr = spl[0]
    low = int(spl[1])
    high = int(spl[2])
    if res.has_key(chr):
      res[chr].insert(low, high)
    else:
      res[chr] = IntervalSet()
      res[chr].insert(low, high)
  return res
  
def intersect_pos(ref, pos):
  res = []
  for chr in sorted(pos.keys()):
    for p in pos[chr]:
      try:
        interval = ref[chr].search(p)
      except KeyError:
        interval = None
      if interval:
        res.append((chr,) + interval)
  return set(res)



if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("reference_bed", type=argparse.FileType("r"))
  parser.add_argument("vcf", 
                      type=argparse.FileType("r"),
                      nargs="+")
  args = parser.parse_args()
 
  bed_intset = make_bed_intset(args.reference_bed)
  args.reference_bed.seek(0)
  all_pos = []
  for l in args.reference_bed:
    spl = l.split()
    all_pos.append((spl[0], int(spl[1]), int(spl[2])))
  
  print "id" + "\t" + "\t".join("%s_%i_%i" % t for t in all_pos)
  for vcf in args.vcf:
    sample_name = vcf.name
    positions = read_vcf(vcf)
    int_pos = intersect_pos(bed_intset, positions)
    flags = [ "2" if p in int_pos else "1" for p in all_pos ]
    print sample_name + "\t" + "\t".join(flags)

