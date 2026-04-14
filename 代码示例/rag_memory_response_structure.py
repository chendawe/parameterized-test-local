# 异步流式（需要 async/await）
async def async_stream_qa(query):
    print("流式响应前的历史：", memory.load_memory_variables({}))
    response_chunks = []
    # 异步 for 循环，必须加 async
    async for chunk in qa_chain.astream({"query": query}):
        if "result" in chunk:
            response_chunks.append(chunk["result"])
            print("当前流式片段：", chunk["result"])
    print("\n流式响应后的历史：", memory.load_memory_variables({}))
    return "".join(response_chunks)

# 调用异步函数
import asyncio
asyncio.run(async_stream_qa("Neo4j 基于什么语言开发？"))


# 自定义记忆键名
memory = ConversationBufferMemory(
    memory_key="my_history",  # 改为自定义键名
    return_messages=True
)
# 此时返回结果中记忆信息的键名就是 "my_history"，而非默认的 "chat_history"
# 典型返回结构（非流式）
{
    "query": "Neo4j 基于什么语言开发？",  # 用户当前提问
    "result": "Neo4j 的核心基于 Java 开发",  # 模型最终回答
    "source_documents": [  # RAG 检索到的原始文档（核心 RAG 信息）
        Document(
            page_content="Neo4j 核心基于 Java 开发",  # 文档内容
            metadata={"source": "xxx", "page": 1}  # 文档元数据（来源、页码等）
        ),
        Document(
            page_content="Neo4j 支持 Bolt 协议",
            metadata={"source": "xxx", "page": 2}
        )
    ],
    "chat_history": [  # 记忆信息（由传入的 memory 提供）
        HumanMessage(content="Neo4j 基于什么语言开发？"),  # 用户消息
        AIMessage(content="Neo4j 的核心基于 Java 开发")  # 模型消息
    ]
}




# 查看prompt模板：
from langchain.chains import RetrievalQA
from langchain.chat_models import ChatOpenAI
from langchain.vectorstores import FAISS
from langchain.embeddings import OpenAIEmbeddings

# 初始化基础组件
llm = ChatOpenAI(model="gpt-3.5-turbo")
embeddings = OpenAIEmbeddings()
vector_db = FAISS.from_texts(["Neo4j 核心基于 Java 开发"], embeddings)
retriever = vector_db.as_retriever()

# 初始化 RetrievalQA 并查看内置 Prompt
qa_chain = RetrievalQA.from_chain_type(
    llm=llm,
    chain_type="stuff",
    retriever=retriever
)

# 打印内置的 Prompt 模板
print("内置 Prompt 模板：")
print(qa_chain.combine_documents_chain.llm_chain.prompt.template)
# 使用以下上下文来回答最后的问题。如果你不知道答案，就说你不知道，不要编造答案。

# {context}  # 这里会填充 RAG 检索到的 source_documents 内容
# {chat_history}  # 这里会填充记忆中的历史对话（如果传入了 memory）

# 问题：{question}
# 答案：


# 自定义模板：
from langchain.prompts import PromptTemplate

# 自定义 Prompt 模板
custom_prompt = PromptTemplate(
    template="""
    你是一个专业的图数据库顾问，必须基于以下上下文和历史对话回答问题，回答要简洁、准确：
    
    历史对话：
    {chat_history}
    
    参考文档：
    {context}
    
    用户问题：{question}
    
    你的回答：
    """,
    input_variables=["chat_history", "context", "question"]
)

# 传入自定义 Prompt 初始化 RetrievalQA
qa_chain = RetrievalQA.from_chain_type(
    llm=llm,
    chain_type="stuff",
    retriever=retriever,
    chain_type_kwargs={"prompt": custom_prompt}  # 指定自定义 Prompt
)






# 完整示例：RAG + 持久化记忆（带 thread_id）
# 以 Redis 为例（需先安装 Redis 和依赖：pip install redis langchain-redis），展示 thread_id 的作用：

from langchain.chains import RetrievalQA
from langchain.chat_models import ChatOpenAI
from langchain.memory import ConversationBufferMemory
from langchain_redis import RedisChatMessageHistory
from langchain.vectorstores import FAISS
from langchain.embeddings import OpenAIEmbeddings
import os

os.environ["OPENAI_API_KEY"] = "your-api-key"
os.environ["REDIS_URL"] = "redis://localhost:6379/0"  # Redis 连接地址

# 1. 初始化向量库（RAG 检索源）
embeddings = OpenAIEmbeddings()
vector_db = FAISS.from_texts(["Neo4j 核心基于 Java 开发"], embeddings)
retriever = vector_db.as_retriever()

# 2. 为不同会话绑定 thread_id（核心：区分用户/会话）
def get_memory(thread_id: str):
    # 初始化 Redis 持久化的会话历史（通过 thread_id 隔离）
    message_history = RedisChatMessageHistory(
        session_id=thread_id,  # 对应你说的 thread_id，唯一标识会话
        redis_url=os.environ["REDIS_URL"]
    )
    # 把持久化的历史传入记忆组件
    memory = ConversationBufferMemory(
        memory_key="chat_history",
        chat_memory=message_history,  # 绑定持久化的历史
        return_messages=True
    )
    return memory

# 3. 初始化 RetrievalQA（针对不同会话）
def get_qa_chain(thread_id: str):
    llm = ChatOpenAI(model="gpt-3.5-turbo", streaming=True)
    qa_chain = RetrievalQA.from_chain_type(
        llm=llm,
        chain_type="stuff",
        retriever=retriever,
        memory=get_memory(thread_id),  # 传入带 thread_id 的记忆
        return_source_documents=True
    )
    return qa_chain

# 4. 测试多会话隔离（不同 thread_id 对应不同记忆）
# 会话 1：thread_id = "user_001"
qa_chain_1 = get_qa_chain("user_001")
qa_chain_1({"query": "Neo4j 基于什么语言开发？"})

# 会话 2：thread_id = "user_002"（不会看到 user_001 的历史）
qa_chain_2 = get_qa_chain("user_002")
qa_chain_2({"query": "Neo4j 支持什么协议？"})

# 重启程序后，通过 thread_id 可恢复对应会话的记忆
qa_chain_1_restore = get_qa_chain("user_001")
print("恢复 user_001 的历史：", qa_chain_1_restore.memory.load_memory_variables({}))








# 持久化记忆
from langchain.chains import RetrievalQA
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain.vectorstores import FAISS
from langchain.memory import ConversationBufferMemory
from langchain_redis import RedisChatMessageHistory
from langchain_core.runnables import RunnableConfig
import os

os.environ["OPENAI_API_KEY"] = "你的API密钥"
os.environ["REDIS_URL"] = "redis://localhost:6379/0"

# 1. 初始化 RAG 检索组件
embeddings = OpenAIEmbeddings()
vector_db = FAISS.from_texts(["Neo4j 核心基于 Java 开发"], embeddings)
retriever = vector_db.as_retriever()

# 2. 动态记忆工厂：读取自定义的 session_id（替代 thread_id）
def get_dynamic_memory(config: RunnableConfig):
    # 这里读取的是 session_id，而非 thread_id
    session_id = config["configurable"]["session_id"]
    # 用 session_id 作为 Redis 的会话标识
    message_history = RedisChatMessageHistory(
        session_id=session_id,  # 核心：和读取的参数名一致即可
        redis_url=os.environ["REDIS_URL"]
    )
    return ConversationBufferMemory(
        memory_key="chat_history",
        chat_memory=message_history,
        return_messages=True
    )

# 3. 初始化链
llm = ChatOpenAI(model="gpt-3.5-turbo")
qa_chain = RetrievalQA.from_chain_type(
    llm=llm,
    chain_type="stuff",
    retriever=retriever,
    return_source_documents=True
)

# 4. 调用时传入自定义的 session_id（而非 thread_id）
def rag_with_custom_id(query: str, session_id: str):
    config = RunnableConfig(
        configurable={"session_id": session_id}  # 自定义参数名：session_id
    )
    qa_chain.memory = get_dynamic_memory(config)
    return qa_chain.invoke({"query": query}, config=config)

# 测试：用 session_id = "user_001" 调用
result = rag_with_custom_id("Neo4j 基于什么语言开发？", "user_001")
print(result["result"])