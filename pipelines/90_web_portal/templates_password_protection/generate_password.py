f = open("fetch_category.php", "r")
content = f.readlines()
f.close()
import string
import random
material = string.ascii_lowercase+string.ascii_uppercase+"123456789"
pswd = "".join([random.choice(material) for k in range(20)])
with open("password.txt", "w") as pswd_file: pswd_file.write(pswd);
pswd = '\t$password = "{pswd}";\n'.format(pswd = pswd)
content[2] = pswd
content = "".join(content)
with open("fetch_category.php", "w") as new_fetch:new_fetch.write(content);
