import sys

class Config(object):
  """ Class to hold all state information about a run, 
  can be used to restart after a break point 
  """
  def __init__(self,infile):
    """ inits the object and reads a file with state info """
    self._file = infile
    self._categories,self._d,self._key_lists = self._read_config_file()
    
  def _read_config_file(self):
    """ read key-value pairs from the state file self._file """
    d = {}
    key_lists = {}
    categories = []
    category = "NONE"
    with open(self._file) as cfg:
       for line in cfg:
          line = line.strip("\n")
          line = line.strip()
          if not (line.startswith("#") or line == ""):
            if line.startswith("[") and line.endswith("]"):
              category = line.strip("[").strip("]").strip()
              categories.append(category)
              d[category] = {}
              key_lists[category] = []
            else:
              words = line.split("=")
              words = [x.strip() for x in words]
              if len(words) < 2:
                continue;
              temp_d = d.get(category,{})
              temp_d[words[0]] = words[1]
              d[category] = temp_d
              temp_l = key_lists.get(category,[])
              temp_l.append(words[0])
              key_lists[category] = temp_l
    return categories,d,key_lists
 
  def _write_config_file(self):
    """ write key-value pairs from self._d to the state file self.file """
    with open(self._file, "w") as cfg:
      for c in self._categories:
        cfg.write("[" + c + "]")
        cfg.write("\n")
        for k in self._key_lists[c]:
          cfg.write(k)
          cfg.write("=")
          cfg.write(self._d[c][k])
          cfg.write("\n")
  
  def get(self,category,key):
    """ 
    get an dictionary entry named "key" from self._d
    """
    return self._d[category][key]

  def set(self,category,key,value):
    """ 
    set an dictionary entry named "key" from self.d to "value"
    syncronize with file on HDD
    """
    if not category in self._d:
      self._d[category] = {}
      self._key_lists[category] = []
      self._categories.append(category)
    if not key in self._d[category]:
      self._key_lists[category].append(key)
    self._d[category][key] = value
    self._write_config_file()
  
  def get_categories(self):
    return self._categories

  def get_keys(self,category):
    return self._key_lists[category]

def main():
  test_obj = Config(sys.argv[1])
  for c in test_obj.get_categories():
    print("Category:", c)
    for i in test_obj.get_keys(c):
      print(i, "=", test_obj.get(c,i))
  test_obj.set("abc","zzz","1")
  for c in test_obj.get_categories():
    print("Category:", c)
    for i in test_obj.get_keys(c):
      print(i, "=", test_obj.get(c, i))


if __name__ == "__main__":
  main()

