
def text_create(name, msg):
    desktop_path = '/Txt/'
    full_path = desktop_path + name + '.Txt'  # 也可以创建一个.doc的word文档
    file = open(full_path, 'w')
    file.write(msg)
    file.write('\r\n')
